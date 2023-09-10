module CNVegStateType

  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
  use decompMod      , only : bounds_type
  use abortutils     , only : endrun
  use spmdMod        , only : masterproc
  use clm_varpar     , only : nlevsno, nlevgrnd, nlevlak, nlevsoi
  use clm_varctl     , only : use_cn, iulog, fsurdat, use_crop, use_cndv, use_crop_agsys
  use clm_varcon     , only : spval, ispval, grlnd
  use landunit_varcon, only : istsoil, istcrop
  use LandunitType   , only : lun
  use ColumnType     , only : col
  use PatchType      , only : patch
  use AnnualFluxDribbler, only : annual_flux_dribbler_type, annual_flux_dribbler_patch
  use dynSubgridControlMod, only : get_for_testing_allow_non_annual_changes
  use CropReprPoolsMod, only : nrepr
  !
  ! !PUBLIC TYPES:
  implicit none
  private
  !
  ! !PUBLIC TYPES:
  type, public :: cnveg_state_type

     integer  , pointer :: burndate_patch              (:)     ! patch crop burn date
     type(annual_flux_dribbler_type) :: dwt_dribbler_patch     ! object to convert instantaneous dwt values into values that are smoothed / dribbled throughout the year
     real(r8) , pointer :: dwt_smoothed_patch          (:)     ! change in patch weight (-1 to 1) on the gridcell in this time step; changes in first time step of year are smoothed (dribbled) over the whole year

     ! Prognostic crop model
     !
     ! TODO(wjs, 2016-02-22) Most / all of these crop-specific state variables should be
     ! moved to CropType
     real(r8) , pointer :: hdidx_patch                 (:)     ! patch cold hardening index?
     real(r8) , pointer :: cumvd_patch                 (:)     ! patch cumulative vernalization d?ependence?
     real(r8) , pointer :: gddmaturity_patch           (:)     ! patch growing degree days (gdd) needed to harvest (ddays)
     real(r8) , pointer :: gddmaturity_thisyr          (:,:)   ! all at-harvest values of the above for this patch this year (ddays) [patch, mxharvests]
     real(r8) , pointer :: huileaf_patch               (:)     ! patch heat unit index needed from planting to leaf emergence
     real(r8) , pointer :: huigrain_patch              (:)     ! patch heat unit index needed to reach vegetative maturity
     real(r8) , pointer :: aleafi_patch                (:)     ! patch saved leaf allocation coefficient from phase 2
     real(r8) , pointer :: astemi_patch                (:)     ! patch saved stem allocation coefficient from phase 2

     real(r8) , pointer :: aleaf_patch                 (:)     ! patch leaf allocation coefficient
     real(r8) , pointer :: astem_patch                 (:)     ! patch stem allocation coefficient
     real(r8) , pointer :: aroot_patch                 (:)     ! patch root allocation coefficient
     real(r8) , pointer :: arepr_patch                 (:,:)   ! patch reproductive allocation coefficient(s)

     ! The following nitrogen-based allocation fractions are just used with the AgSys
     ! crop model (use_crop_agsys = .true.):
     real(r8) , pointer :: aleaf_n_patch               (:)     ! patch leaf allocation coefficient for N
     real(r8) , pointer :: astem_n_patch               (:)     ! patch stem allocation coefficient for N
     real(r8) , pointer :: aroot_n_patch               (:)     ! patch root allocation coefficient for N
     real(r8) , pointer :: arepr_n_patch               (:,:)   ! patch reproductive allocation coefficient(s) for N

     real(r8) , pointer :: htmx_patch                  (:)     ! patch max hgt attained by a crop during yr (m)
     integer  , pointer :: peaklai_patch               (:)     ! patch 1: max allowed lai; 0: not at max

     integer  , pointer :: idop_patch                  (:)     ! patch date of planting (day of year)
     integer  , pointer :: iyop_patch                  (:)     ! patch year of planting

     real(r8) , pointer :: lgdp_col                    (:)     ! col gdp limitation factor for fire occurrence (0-1)
     real(r8) , pointer :: lgdp1_col                   (:)     ! col gdp limitation factor for fire spreading (0-1)
     real(r8) , pointer :: lpop_col                    (:)     ! col pop limitation factor for fire spreading (0-1)

     real(r8) , pointer :: tempavg_t2m_patch           (:)     ! patch temporary average 2m air temperature (K)
     real(r8) , pointer :: annavg_t2m_patch            (:)     ! patch annual average 2m air temperature (K)
     real(r8) , pointer :: annavg_t2m_col              (:)     ! col annual average of 2m air temperature, averaged from patch-level (K)
     real(r8) , pointer :: annsum_counter_col          (:)     ! col seconds since last annual accumulator turnover

     ! Fire
     real(r8) , pointer :: nfire_col                   (:)     ! col fire counts (count/km2/sec), valid only in Reg. C
     real(r8) , pointer :: fsr_col                     (:)     ! col fire spread rate at column level (m/s)
     real(r8) , pointer :: fd_col                      (:)     ! col fire duration at column level (hr)
     real(r8) , pointer :: lfc_col                     (:)     ! col conversion area fraction of BET and BDT that haven't burned before (/timestep)
     real(r8) , pointer :: lfc2_col                    (:)     ! col conversion area fraction of BET and BDT that burned (/sec)
     real(r8) , pointer :: dtrotr_col                  (:)     ! col annual decreased fraction coverage of BET on the gridcell (0-1)
     real(r8) , pointer :: trotr1_col                  (:)     ! col patch weight of BET on the column (0-1)
     real(r8) , pointer :: trotr2_col                  (:)     ! col patch weight of BDT on the column (0-1)
     real(r8) , pointer :: cropf_col                   (:)     ! col crop fraction in veg column (0-1)
     real(r8) , pointer :: baf_crop_col                (:)     ! col baf for cropland(/sec)
     real(r8) , pointer :: baf_peatf_col               (:)     ! col baf for peatland (/sec)
     real(r8) , pointer :: fbac_col                    (:)     ! col total burned area out of conversion (/sec)
     real(r8) , pointer :: fbac1_col                   (:)     ! col burned area out of conversion region due to land use fire (/sec)
     real(r8) , pointer :: wtlf_col                    (:)     ! col fractional coverage of non-crop Patches (0-1)
     real(r8) , pointer :: lfwt_col                    (:)     ! col fractional coverage of non-crop and non-bare-soil Patches (0-1)
     real(r8) , pointer :: farea_burned_col            (:)     ! col fractional area burned (/sec)

     real(r8), pointer :: dormant_flag_patch           (:)     ! patch dormancy flag
     real(r8), pointer :: days_active_patch            (:)     ! patch number of days since last dormancy
     real(r8), pointer :: onset_flag_patch             (:)     ! patch onset flag
     real(r8), pointer :: onset_counter_patch          (:)     ! patch onset days counter
     real(r8), pointer :: onset_gddflag_patch          (:)     ! patch onset flag for growing degree day sum
     real(r8), pointer :: onset_fdd_patch              (:)     ! patch onset freezing degree days counter
     real(r8), pointer :: onset_gdd_patch              (:)     ! patch onset growing degree days
     real(r8), pointer :: onset_swi_patch              (:)     ! patch onset soil water index
     real(r8), pointer :: offset_flag_patch            (:)     ! patch offset flag
     real(r8), pointer :: offset_counter_patch         (:)     ! patch offset days counter
     real(r8), pointer :: offset_fdd_patch             (:)     ! patch offset freezing degree days counter
     real(r8), pointer :: offset_swi_patch             (:)     ! patch offset soil water index
     real(r8), pointer :: grain_flag_patch             (:)     ! patch 1: grain fill stage; 0: not
     real(r8), pointer :: lgsf_patch                   (:)     ! patch long growing season factor [0-1]
     real(r8), pointer :: bglfr_patch                  (:)     ! patch background litterfall rate (1/s)
     real(r8), pointer :: bgtr_patch                   (:)     ! patch background transfer growth rate (1/s)
     real(r8), pointer :: c_allometry_patch            (:)     ! patch C allocation index (DIM)
     real(r8), pointer :: n_allometry_patch            (:)     ! patch N allocation index (DIM)

     real(r8), pointer :: tempsum_potential_gpp_patch  (:)     ! patch temporary annual sum of potential GPP
     real(r8), pointer :: annsum_potential_gpp_patch   (:)     ! patch annual sum of potential GPP
     real(r8), pointer :: tempmax_retransn_patch       (:)     ! patch temporary annual max of retranslocated N pool (gN/m2)
     real(r8), pointer :: annmax_retransn_patch        (:)     ! patch annual max of retranslocated N pool (gN/m2)
     real(r8), pointer :: downreg_patch                (:)     ! patch fractional reduction in GPP due to N limitation (DIM)
     real(r8), pointer :: leafcn_offset_patch          (:)     ! patch leaf C:N used by FUN
     real(r8), pointer :: plantCN_patch                (:)     ! patch plant C:N used by FUN

   contains

     procedure, public  :: Init
     procedure, public  :: Restart
     procedure, private :: InitAllocate
     procedure, private :: InitHistory
     procedure, private :: InitCold

  end type cnveg_state_type
  !------------------------------------------------------------------------

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds, alloc_full_veg)

    class(cnveg_state_type) :: this
    type(bounds_type), intent(in) :: bounds
    logical,intent(in) :: alloc_full_veg  ! Total number of bgc patches on proc (non-fates)

    call this%InitAllocate ( bounds, alloc_full_veg)
    if (use_cn) then
       call this%InitHistory ( bounds )
    end if
    if(alloc_full_veg) then  !This is true if not use_fates_bgc
       call this%InitCold ( bounds )
    end if
    
  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds, alloc_full_veg)
    !
    ! !DESCRIPTION:
    ! Initialize module data structure
    !
    ! !USES:
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
    use clm_varpar, only : mxsowings, mxharvests
    !
    ! !ARGUMENTS:
    class(cnveg_state_type) :: this
    type(bounds_type), intent(in) :: bounds
    logical, intent(in) :: alloc_full_veg ! Total number of bgc patches on proc (non-fates)
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    integer :: begc, endc
    logical :: allows_non_annual_delta
    !------------------------------------------------------------------------

    if(alloc_full_veg)then
       begp = bounds%begp; endp= bounds%endp
       begc = bounds%begc; endc= bounds%endc
    else
       begp = 0;endp = 0
       begc = 0;endc = 0
    end if
       
       
    ! Note that we set allows_non_annual_delta to false because we expect land cover
    ! change to be applied entirely at the start of the year. Currently the fire code
    ! appears to assume that the land cover change rate is constant throughout the year,
    ! in this code (which is accompanied by the comment, 'land cover conversion in CLM4.5
    ! is the same for each timestep except for the beginning'):
    !
    !         if( kmo == 1 .and. kda == 1 .and. mcsec == dt)then
    !            lfc(c) = dtrotr_col(c)*dayspyr*secspday/dt
    !         end if
    !
    ! so setting allows_non_annual_delta to .false. helps ensure that remains true.
    !
    ! However, we do keep allows_non_annual_delta = .true. if running with CNDV, because
    ! (in contrast with other land cover change) CNDV currently still interpolates land
    ! cover change throughout the year. Note that there is therefore an inconsistency with
    ! the fire code if we're using CNDV, due to the way the annual flux dribbler works:
    ! The dwt generated by CNDV on the first time step of the year is dribbled throughout
    ! the year by dwt_dribbler_patch, but the CNDV dwt on every other time step comes at
    ! its full value. So there will be a lower dwt in the first time step of the year
    ! relative to every other time step of the year. If CNDV is the main contributor to
    ! dwt, I think this can lead to a large violation of the above assumption of constant
    ! dwt in the fire code. However, the fire code doesn't seem designed to work with
    ! CNDV at all (because land cover change is assumed to be associated with
    ! deforestation, not natural changes in areas), so maybe this inconsistency is the
    ! least of the problem: see bug 2392.
    if (get_for_testing_allow_non_annual_changes()) then
       allows_non_annual_delta = .true.
    else if (use_cndv) then
       allows_non_annual_delta = .true.
    else
       allows_non_annual_delta = .false.
    end if
    this%dwt_dribbler_patch = annual_flux_dribbler_patch( &
         bounds = bounds, &
         name = 'dwt', &
         units = 'fractional area', &
         allows_non_annual_delta = allows_non_annual_delta)

    allocate(this%burndate_patch      (begp:endp))                   ; this%burndate_patch      (:)   = ispval
    allocate(this%dwt_smoothed_patch  (begp:endp))                   ; this%dwt_smoothed_patch  (:)   = nan

    allocate(this%hdidx_patch         (begp:endp))                   ; this%hdidx_patch         (:)   = nan
    allocate(this%cumvd_patch         (begp:endp))                   ; this%cumvd_patch         (:)   = nan
    allocate(this%gddmaturity_patch   (begp:endp))                   ; this%gddmaturity_patch   (:)   = spval
    allocate(this%gddmaturity_thisyr  (begp:endp,1:mxharvests))      ; this%gddmaturity_thisyr  (:,:) = spval
    allocate(this%huileaf_patch       (begp:endp))                   ; this%huileaf_patch       (:)   = nan
    allocate(this%huigrain_patch      (begp:endp))                   ; this%huigrain_patch      (:)   = 0.0_r8
    allocate(this%aleafi_patch        (begp:endp))                   ; this%aleafi_patch        (:)   = nan
    allocate(this%astemi_patch        (begp:endp))                   ; this%astemi_patch        (:)   = nan

    allocate(this%aleaf_patch         (begp:endp))                   ; this%aleaf_patch         (:)   = nan
    allocate(this%astem_patch         (begp:endp))                   ; this%astem_patch         (:)   = nan
    allocate(this%aroot_patch         (begp:endp))                   ; this%aroot_patch         (:)   = nan
    allocate(this%arepr_patch         (begp:endp, nrepr))            ; this%arepr_patch         (:,:) = nan

    if (use_crop_agsys) then
       allocate(this%aleaf_n_patch    (begp:endp))                   ; this%aleaf_n_patch       (:)   = nan
       allocate(this%astem_n_patch    (begp:endp))                   ; this%astem_n_patch       (:)   = nan
       allocate(this%aroot_n_patch    (begp:endp))                   ; this%aroot_n_patch       (:)   = nan
       allocate(this%arepr_n_patch    (begp:endp, nrepr))            ; this%arepr_n_patch       (:,:) = nan
    end if

    allocate(this%htmx_patch          (begp:endp))                   ; this%htmx_patch          (:)   = 0.0_r8
    allocate(this%peaklai_patch       (begp:endp))                   ; this%peaklai_patch       (:)   = 0

    allocate(this%idop_patch          (begp:endp))                   ; this%idop_patch          (:)   = huge(1)
    allocate(this%iyop_patch          (begp:endp))                   ; this%iyop_patch          (:)   = ispval

    allocate(this%lgdp_col            (begc:endc))                   ;
    allocate(this%lgdp1_col           (begc:endc))                   ;
    allocate(this%lpop_col            (begc:endc))                   ;

    allocate(this%tempavg_t2m_patch   (begp:endp))                   ; this%tempavg_t2m_patch   (:)   = nan
    allocate(this%annsum_counter_col  (begc:endc))                   ; this%annsum_counter_col  (:)   = nan
    allocate(this%annavg_t2m_col      (begc:endc))                   ; this%annavg_t2m_col      (:)   = nan
    allocate(this%annavg_t2m_patch    (begp:endp))                   ; this%annavg_t2m_patch    (:)   = nan

    allocate(this%nfire_col           (begc:endc))                   ; this%nfire_col           (:)   = spval
    allocate(this%fsr_col             (begc:endc))                   ; this%fsr_col             (:)   = nan
    allocate(this%fd_col              (begc:endc))                   ; this%fd_col              (:)   = nan
    allocate(this%lfc_col             (begc:endc))                   ; this%lfc_col             (:)   = spval
    allocate(this%lfc2_col            (begc:endc))                   ; this%lfc2_col            (:)   = 0._r8
    allocate(this%dtrotr_col          (begc:endc))                   ; this%dtrotr_col          (:)   = 0._r8
    allocate(this%trotr1_col          (begc:endc))                   ; this%trotr1_col          (:)   = 0._r8
    allocate(this%trotr2_col          (begc:endc))                   ; this%trotr2_col          (:)   = 0._r8
    allocate(this%cropf_col           (begc:endc))                   ; this%cropf_col           (:)   = nan
    allocate(this%baf_crop_col        (begc:endc))                   ; this%baf_crop_col        (:)   = nan
    allocate(this%baf_peatf_col       (begc:endc))                   ; this%baf_peatf_col       (:)   = nan
    allocate(this%fbac_col            (begc:endc))                   ; this%fbac_col            (:)   = nan
    allocate(this%fbac1_col           (begc:endc))                   ; this%fbac1_col           (:)   = nan
    allocate(this%wtlf_col            (begc:endc))                   ; this%wtlf_col            (:)   = nan
    allocate(this%lfwt_col            (begc:endc))                   ; this%lfwt_col            (:)   = nan
    allocate(this%farea_burned_col    (begc:endc))                   ; this%farea_burned_col    (:)   = nan

    allocate(this%dormant_flag_patch          (begp:endp)) ;    this%dormant_flag_patch          (:) = nan
    allocate(this%days_active_patch           (begp:endp)) ;    this%days_active_patch           (:) = nan
    allocate(this%onset_flag_patch            (begp:endp)) ;    this%onset_flag_patch            (:) = nan
    allocate(this%onset_counter_patch         (begp:endp)) ;    this%onset_counter_patch         (:) = nan
    allocate(this%onset_gddflag_patch         (begp:endp)) ;    this%onset_gddflag_patch         (:) = nan
    allocate(this%onset_fdd_patch             (begp:endp)) ;    this%onset_fdd_patch             (:) = nan
    allocate(this%onset_gdd_patch             (begp:endp)) ;    this%onset_gdd_patch             (:) = nan
    allocate(this%onset_swi_patch             (begp:endp)) ;    this%onset_swi_patch             (:) = nan
    allocate(this%offset_flag_patch           (begp:endp)) ;    this%offset_flag_patch           (:) = nan
    allocate(this%offset_counter_patch        (begp:endp)) ;    this%offset_counter_patch        (:) = nan
    allocate(this%offset_fdd_patch            (begp:endp)) ;    this%offset_fdd_patch            (:) = nan
    allocate(this%offset_swi_patch            (begp:endp)) ;    this%offset_swi_patch            (:) = nan
    allocate(this%grain_flag_patch            (begp:endp)) ;    this%grain_flag_patch            (:) = nan
    allocate(this%lgsf_patch                  (begp:endp)) ;    this%lgsf_patch                  (:) = nan
    allocate(this%bglfr_patch                 (begp:endp)) ;    this%bglfr_patch                 (:) = nan
    allocate(this%bgtr_patch                  (begp:endp)) ;    this%bgtr_patch                  (:) = nan
    allocate(this%c_allometry_patch           (begp:endp)) ;    this%c_allometry_patch           (:) = nan
    allocate(this%n_allometry_patch           (begp:endp)) ;    this%n_allometry_patch           (:) = nan
    allocate(this%tempsum_potential_gpp_patch (begp:endp)) ;    this%tempsum_potential_gpp_patch (:) = nan
    allocate(this%annsum_potential_gpp_patch  (begp:endp)) ;    this%annsum_potential_gpp_patch  (:) = nan
    allocate(this%tempmax_retransn_patch      (begp:endp)) ;    this%tempmax_retransn_patch      (:) = nan
    allocate(this%annmax_retransn_patch       (begp:endp)) ;    this%annmax_retransn_patch       (:) = nan
    allocate(this%downreg_patch               (begp:endp)) ;    this%downreg_patch               (:) = nan
    allocate(this%leafcn_offset_patch         (begp:endp)) ;    this%leafcn_offset_patch         (:) = nan
    allocate(this%plantCN_patch               (begp:endp)) ;    this%plantCN_patch               (:) = nan

  end subroutine InitAllocate

  !------------------------------------------------------------------------
  subroutine InitHistory(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize module data structure
    !
    ! !USES:
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
    use histFileMod    , only : hist_addfld1d, hist_addfld2d, hist_addfld_decomp, no_snow_normal
    !
    ! !ARGUMENTS:
    class(cnveg_state_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer           :: begp, endp
    integer           :: begc, endc
    character(8)      :: vr_suffix
    character(10)     :: active
    real(r8), pointer :: data2dptr(:,:), data1dptr(:) ! temp. pointers for slicing larger arrays
    !------------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp
    begc = bounds%begc; endc= bounds%endc

    if ( use_crop) then
       ! Daily
       this%gddmaturity_patch(begp:endp) = spval
       call hist_addfld1d (fname='GDDHARV', units='ddays', &
            avgflag='A', long_name='Growing degree days (gdd) needed to harvest', &
            ptr_patch=this%gddmaturity_patch, default='inactive')
       
       ! Per harvest
       this%gddmaturity_thisyr(begp:endp,:) = spval
       call hist_addfld2d (fname='GDDHARV_PERHARV', units='ddays', type2d='mxharvests', &
            avgflag='I', long_name='Growing degree days (gdd) needed to harvest; should only be output annually', &
            ptr_patch=this%gddmaturity_thisyr, default='inactive')

    end if

    this%lfc2_col(begc:endc) = spval
    call hist_addfld1d (fname='LFC2', units='per sec', &
         avgflag='A', long_name='conversion area fraction of BET and BDT that burned', &
         ptr_col=this%lfc2_col)

    this%annsum_counter_col(begc:endc) = spval
    call hist_addfld1d (fname='ANNSUM_COUNTER', units='s', &
         avgflag='A', long_name='seconds since last annual accumulator turnover', &
         ptr_col=this%annsum_counter_col, default='inactive')

    this%annavg_t2m_col(begc:endc) = spval
    call hist_addfld1d (fname='CANNAVG_T2M', units='K', &
         avgflag='A', long_name='annual average of 2m air temperature', &
         ptr_col=this%annavg_t2m_col, default='inactive')

    this%nfire_col(begc:endc) = spval
    call hist_addfld1d (fname='NFIRE',  units='counts/km2/sec', &
         avgflag='A', long_name='fire counts valid only in Reg.C', &
         ptr_col=this%nfire_col)

    this%farea_burned_col(begc:endc) = spval
    call hist_addfld1d (fname='FAREA_BURNED',  units='s-1', &
         avgflag='A', long_name='timestep fractional area burned', &
         ptr_col=this%farea_burned_col)

    this%baf_crop_col(begc:endc) = spval
    call hist_addfld1d (fname='BAF_CROP',  units='s-1', &
         avgflag='A', long_name='fractional area burned for crop', &
         ptr_col=this%baf_crop_col)

    this%baf_peatf_col(begc:endc) = spval
    call hist_addfld1d (fname='BAF_PEATF',  units='s-1', &
         avgflag='A', long_name='fractional area burned in peatland', &
         ptr_col=this%baf_peatf_col)

    this%annavg_t2m_patch(begp:endp) = spval
    call hist_addfld1d (fname='ANNAVG_T2M', units='K', &
         avgflag='A', long_name='annual average 2m air temperature', &
         ptr_patch=this%annavg_t2m_patch, default='inactive')

    this%tempavg_t2m_patch(begp:endp) = spval
    call hist_addfld1d (fname='TEMPAVG_T2M', units='K', &
         avgflag='A', long_name='temporary average 2m air temperature', &
         ptr_patch=this%tempavg_t2m_patch, default='inactive')

    this%dormant_flag_patch(begp:endp) = spval
    call hist_addfld1d (fname='DORMANT_FLAG', units='none', &
         avgflag='A', long_name='dormancy flag', &
         ptr_patch=this%dormant_flag_patch, default='inactive')

    this%days_active_patch(begp:endp) = spval
    call hist_addfld1d (fname='DAYS_ACTIVE', units='days', &
         avgflag='A', long_name='number of days since last dormancy', &
         ptr_patch=this%days_active_patch, default='inactive')

    this%onset_flag_patch(begp:endp) = spval
    call hist_addfld1d (fname='ONSET_FLAG', units='none', &
         avgflag='A', long_name='onset flag', &
         ptr_patch=this%onset_flag_patch, default='inactive')

    this%onset_counter_patch(begp:endp) = spval
    call hist_addfld1d (fname='ONSET_COUNTER', units='days', &
         avgflag='A', long_name='onset days counter', &
         ptr_patch=this%onset_counter_patch, default='inactive')

    this%onset_gddflag_patch(begp:endp) = spval
    call hist_addfld1d (fname='ONSET_GDDFLAG', units='none', &
         avgflag='A', long_name='onset flag for growing degree day sum', &
         ptr_patch=this%onset_gddflag_patch, default='inactive')

    this%onset_fdd_patch(begp:endp) = spval
    call hist_addfld1d (fname='ONSET_FDD', units='C degree-days', &
         avgflag='A', long_name='onset freezing degree days counter', &
         ptr_patch=this%onset_fdd_patch, default='inactive')

    this%onset_gdd_patch(begp:endp) = spval
    call hist_addfld1d (fname='ONSET_GDD', units='C degree-days', &
         avgflag='A', long_name='onset growing degree days', &
         ptr_patch=this%onset_gdd_patch, default='inactive')

    this%onset_swi_patch(begp:endp) = spval
    call hist_addfld1d (fname='ONSET_SWI', units='none', &
         avgflag='A', long_name='onset soil water index', &
         ptr_patch=this%onset_swi_patch, default='inactive')

    this%offset_flag_patch(begp:endp) = spval
    call hist_addfld1d (fname='OFFSET_FLAG', units='none', &
         avgflag='A', long_name='offset flag', &
         ptr_patch=this%offset_flag_patch, default='inactive')

    this%offset_counter_patch(begp:endp) = spval
    call hist_addfld1d (fname='OFFSET_COUNTER', units='days', &
         avgflag='A', long_name='offset days counter', &
         ptr_patch=this%offset_counter_patch, default='inactive')

    this%offset_fdd_patch(begp:endp) = spval
    call hist_addfld1d (fname='OFFSET_FDD', units='C degree-days', &
         avgflag='A', long_name='offset freezing degree days counter', &
         ptr_patch=this%offset_fdd_patch, default='inactive')

    this%offset_swi_patch(begp:endp) = spval
    call hist_addfld1d (fname='OFFSET_SWI', units='none', &
         avgflag='A', long_name='offset soil water index', &
         ptr_patch=this%offset_swi_patch, default='inactive')

    this%lgsf_patch(begp:endp) = spval
    call hist_addfld1d (fname='LGSF', units='proportion', &
         avgflag='A', long_name='long growing season factor', &
         ptr_patch=this%lgsf_patch, default='inactive')

    this%bglfr_patch(begp:endp) = spval
    call hist_addfld1d (fname='BGLFR', units='1/s', &
         avgflag='A', long_name='background litterfall rate', &
         ptr_patch=this%bglfr_patch, default='inactive')

    this%bgtr_patch(begp:endp) = spval
    call hist_addfld1d (fname='BGTR', units='1/s', &
         avgflag='A', long_name='background transfer growth rate', &
         ptr_patch=this%bgtr_patch, default='inactive')

    this%c_allometry_patch(begp:endp) = spval
    call hist_addfld1d (fname='C_ALLOMETRY', units='none', &
         avgflag='A', long_name='C allocation index', &
         ptr_patch=this%c_allometry_patch, default='inactive')

    this%n_allometry_patch(begp:endp) = spval
    call hist_addfld1d (fname='N_ALLOMETRY', units='none', &
         avgflag='A', long_name='N allocation index', &
         ptr_patch=this%n_allometry_patch, default='inactive')

    this%tempsum_potential_gpp_patch(begp:endp) = spval
    call hist_addfld1d (fname='TEMPSUM_POTENTIAL_GPP', units='gC/m^2/yr', &
         avgflag='A', long_name='temporary annual sum of potential GPP', &
         ptr_patch=this%tempsum_potential_gpp_patch, default='inactive')

    this%annsum_potential_gpp_patch(begp:endp) = spval
    call hist_addfld1d (fname='ANNSUM_POTENTIAL_GPP', units='gN/m^2/yr', &
         avgflag='A', long_name='annual sum of potential GPP', &
         ptr_patch=this%annsum_potential_gpp_patch, default='inactive')

    this%tempmax_retransn_patch(begp:endp) = spval
    call hist_addfld1d (fname='TEMPMAX_RETRANSN', units='gN/m^2', &
         avgflag='A', long_name='temporary annual max of retranslocated N pool', &
         ptr_patch=this%tempmax_retransn_patch, default='inactive')

    this%annmax_retransn_patch(begp:endp) = spval
    call hist_addfld1d (fname='ANNMAX_RETRANSN', units='gN/m^2', &
         avgflag='A', long_name='annual max of retranslocated N pool', &
         ptr_patch=this%annmax_retransn_patch, default='inactive')

    this%downreg_patch(begp:endp) = spval
    call hist_addfld1d (fname='DOWNREG', units='proportion', &
         avgflag='A', long_name='fractional reduction in GPP due to N limitation', &
         ptr_patch=this%downreg_patch, default='inactive')

    this%leafcn_offset_patch(begp:endp) = spval
    call hist_addfld1d (fname='LEAFCN_OFFSET', units='unitless', &
         avgflag='A', long_name='Leaf C:N used by FUN', &
         ptr_patch=this%leafcn_offset_patch, default='inactive')

    this%plantCN_patch(begp:endp)       = spval
    call hist_addfld1d (fname='PLANTCN', units='unitless', &
         avgflag='A', long_name='Plant C:N used by FUN', &
         ptr_patch=this%plantCN_patch, default='inactive')
  end subroutine InitHistory

  !-----------------------------------------------------------------------
  subroutine InitCold(this, bounds)
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(cnveg_state_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer               :: g,l,c,p            ! dices
    !-----------------------------------------------------------------------

    ! --------------------------------------------------------------------
    ! Initialize terms needed for dust model
    ! TODO - move these terms to DUSTMod module variables
    ! --------------------------------------------------------------------

    do c = bounds%begc, bounds%endc
       l = col%landunit(c)
       if (lun%ifspecial(l)) then
          this%annsum_counter_col (c) = spval
          this%annavg_t2m_col     (c) = spval
          this%nfire_col          (c) = spval
          this%baf_crop_col       (c) = spval
          this%baf_peatf_col      (c) = spval
          this%fbac_col           (c) = spval
          this%fbac1_col          (c) = spval
          this%farea_burned_col   (c) = spval
       end if

       if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then
          this%annsum_counter_col(c) = 0._r8
          this%annavg_t2m_col(c)     = 280._r8

          ! fire related variables
          this%baf_crop_col(c)       = 0._r8
          this%baf_peatf_col(c)      = 0._r8
          this%fbac_col(c)           = 0._r8
          this%fbac1_col(c)          = 0._r8
          this%farea_burned_col(c)   = 0._r8
          this%nfire_col(c)          = 0._r8
       end if
    end do

    ! ecophysiological and phenology variables

    do p = bounds%begp,bounds%endp
       l = patch%landunit(p)

       if (lun%ifspecial(l)) then
          this%annavg_t2m_patch  (p)          = spval
          this%tempavg_t2m_patch (p)          = spval
          this%dormant_flag_patch(p)          = spval
          this%days_active_patch(p)           = spval
          this%onset_flag_patch(p)            = spval
          this%onset_counter_patch(p)         = spval
          this%onset_gddflag_patch(p)         = spval
          this%onset_fdd_patch(p)             = spval
          this%onset_gdd_patch(p)             = spval
          this%onset_swi_patch(p)             = spval
          this%offset_flag_patch(p)           = spval
          this%offset_counter_patch(p)        = spval
          this%offset_fdd_patch(p)            = spval
          this%offset_swi_patch(p)            = spval
          this%grain_flag_patch(p)            = spval
          this%lgsf_patch(p)                  = spval
          this%bglfr_patch(p)                 = spval
          this%bgtr_patch(p)                  = spval
          this%c_allometry_patch(p)           = spval
          this%n_allometry_patch(p)           = spval
          this%tempsum_potential_gpp_patch(p) = spval
          this%annsum_potential_gpp_patch(p)  = spval
          this%tempmax_retransn_patch(p)      = spval
          this%annmax_retransn_patch(p)       = spval
          this%downreg_patch(p)               = spval
          this%leafcn_offset_patch(p)         = spval
          this%plantCN_patch(p)               = spval
       end if

       if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then
          ! phenology variables
          this%dormant_flag_patch(p)   = 1._r8
          this%days_active_patch(p)    = 0._r8
          this%onset_flag_patch(p)     = 0._r8
          this%onset_counter_patch(p)  = 0._r8
          this%onset_gddflag_patch(p)  = 0._r8
          this%onset_fdd_patch(p)      = 0._r8
          this%onset_gdd_patch(p)      = 0._r8
          this%onset_swi_patch(p)      = 0._r8
          this%offset_flag_patch(p)    = 0._r8
          this%offset_counter_patch(p) = 0._r8
          this%offset_fdd_patch(p)     = 0._r8
          this%offset_swi_patch(p)     = 0._r8
          this%lgsf_patch(p)           = 0._r8
          this%bglfr_patch(p)          = 0._r8
          this%bgtr_patch(p)           = 0._r8
          this%annavg_t2m_patch(p)     = 280._r8
          this%tempavg_t2m_patch(p)    = 0._r8
          this%grain_flag_patch(p)     = 0._r8

          ! non-phenology variables
          this%c_allometry_patch(p)           = 0._r8
          this%n_allometry_patch(p)           = 0._r8
          this%tempsum_potential_gpp_patch(p) = 0._r8
          this%annsum_potential_gpp_patch(p)  = 0._r8
          this%tempmax_retransn_patch(p)      = 0._r8
          this%annmax_retransn_patch(p)       = 0._r8
          this%downreg_patch(p)               = 0._r8
          this%leafcn_offset_patch(p)         = spval
          this%plantCN_patch(p)               = spval
       end if

    end do

    ! fire variables

    do c = bounds%begc,bounds%endc
       this%lfc2_col(c) = 0._r8
    end do

  end subroutine InitCold

  !------------------------------------------------------------------------
  subroutine Restart(this, bounds, ncid, flag, cnveg_carbonstate, &
             cnveg_nitrogenstate, filter_reseed_patch, num_reseed_patch)
    !
    ! !USES:
    use shr_log_mod, only : errMsg => shr_log_errMsg
    use spmdMod    , only : masterproc
    use abortutils , only : endrun
    use CNVegNitrogenStateType, only: cnveg_nitrogenstate_type
    use CNVegCarbonStateType  , only: cnveg_carbonstate_type
    use restUtilMod
    use ncdio_pio
    use pftconMod , only : pftcon
    !
    ! !ARGUMENTS:
    class(cnveg_state_type) :: this
    type(bounds_type), intent(in)    :: bounds
    type(file_desc_t), intent(inout) :: ncid
    character(len=*) , intent(in)    :: flag
    type(cnveg_nitrogenstate_type), intent(in) :: cnveg_nitrogenstate
    type(cnveg_carbonstate_type)  , intent(in) :: cnveg_carbonstate
    integer                       , intent(out), optional :: filter_reseed_patch(:)
    integer                       , intent(out), optional :: num_reseed_patch
    !
    ! !LOCAL VARIABLES:
    integer          :: j,c,i,p ! indices
    logical          :: readvar   ! determine if variable is on initial file
    real(r8), pointer :: ptr2d(:,:) ! temp. pointers for slicing larger arrays
    real(r8), pointer :: ptr1d(:)   ! temp. pointers for slicing larger arrays
    !-----------------------------------------------------------------------

    call this%dwt_dribbler_patch%Restart(bounds, ncid, flag)

    call restartvar(ncid=ncid, flag=flag, varname='dormant_flag', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='dormancy flag', units='unitless', &
         interpinic_flag='interp', readvar=readvar, data=this%dormant_flag_patch)

    call restartvar(ncid=ncid, flag=flag, varname='days_active', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='number of days since last dormancy', units='days' , &
         interpinic_flag='interp', readvar=readvar, data=this%days_active_patch)

    call restartvar(ncid=ncid, flag=flag, varname='onset_flag', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='flag if critical growing degree-day sum is exceeded', units='unitless' , &
         interpinic_flag='interp', readvar=readvar, data=this%onset_flag_patch)

    call restartvar(ncid=ncid, flag=flag, varname='onset_counter', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='onset days counter', units='sec' , &
         interpinic_flag='interp', readvar=readvar, data=this%onset_counter_patch)

    call restartvar(ncid=ncid, flag=flag, varname='onset_gddflag', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='onset flag for growing degree day sum', units='' , &
         interpinic_flag='interp', readvar=readvar, data=this%onset_gddflag_patch)

    call restartvar(ncid=ncid, flag=flag, varname='onset_fdd', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='onset freezing degree days counter', units='days' , &
         interpinic_flag='interp', readvar=readvar, data=this%onset_fdd_patch)

    call restartvar(ncid=ncid, flag=flag, varname='onset_gdd', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='onset growing degree days', units='days' , &
         interpinic_flag='interp', readvar=readvar, data=this%onset_gdd_patch)

    call restartvar(ncid=ncid, flag=flag, varname='onset_swi', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='onset soil water index', units='days' , &
         interpinic_flag='interp', readvar=readvar, data=this%onset_swi_patch)

    call restartvar(ncid=ncid, flag=flag, varname='offset_flag', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='offset flag', units='unitless' , &
         interpinic_flag='interp', readvar=readvar, data=this%offset_flag_patch)

    call restartvar(ncid=ncid, flag=flag, varname='offset_counter', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='offset days counter', units='sec' , &
         interpinic_flag='interp', readvar=readvar, data=this%offset_counter_patch)

    call restartvar(ncid=ncid, flag=flag, varname='offset_fdd', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='offset freezing degree days counter', units='days' , &
         interpinic_flag='interp', readvar=readvar, data=this%offset_fdd_patch)

    call restartvar(ncid=ncid, flag=flag, varname='offset_swi', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%offset_swi_patch)

    call restartvar(ncid=ncid, flag=flag, varname='lgsf', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%lgsf_patch)

    call restartvar(ncid=ncid, flag=flag, varname='bglfr', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%bglfr_patch)

    call restartvar(ncid=ncid, flag=flag, varname='bgtr', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%bgtr_patch)

    call restartvar(ncid=ncid, flag=flag, varname='annavg_t2m', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%annavg_t2m_patch)

    call restartvar(ncid=ncid, flag=flag, varname='tempavg_t2m', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%tempavg_t2m_patch)

    call restartvar(ncid=ncid, flag=flag, varname='c_allometry', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%c_allometry_patch)

    call restartvar(ncid=ncid, flag=flag, varname='n_allometry', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%n_allometry_patch)

    call restartvar(ncid=ncid, flag=flag, varname='tempsum_potential_gpp', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%tempsum_potential_gpp_patch)

    call restartvar(ncid=ncid, flag=flag, varname='annsum_potential_gpp', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%annsum_potential_gpp_patch)

    call restartvar(ncid=ncid, flag=flag, varname='tempmax_retransn', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%tempmax_retransn_patch)

    call restartvar(ncid=ncid, flag=flag, varname='annmax_retransn', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%annmax_retransn_patch)

    call restartvar(ncid=ncid, flag=flag, varname='downreg', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%downreg_patch)

    call restartvar(ncid=ncid, flag=flag, varname='leafcn_offset', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%leafcn_offset_patch)

    call restartvar(ncid=ncid, flag=flag, varname='plantCN', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%plantCN_patch)

    call restartvar(ncid=ncid, flag=flag, varname='annsum_counter', xtype=ncd_double,  &
         dim1name='column', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%annsum_counter_col)

    call restartvar(ncid=ncid, flag=flag, varname='burndate', xtype=ncd_int,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%burndate_patch)

    call restartvar(ncid=ncid, flag=flag, varname='lfc', xtype=ncd_double,  &
         dim1name='column', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%lfc_col)

    call restartvar(ncid=ncid, flag=flag, varname='cannavg_t2m', xtype=ncd_double,  &
         dim1name='column', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%annavg_t2m_col)

    if (use_crop) then

       call restartvar(ncid=ncid, flag=flag,  varname='htmx', xtype=ncd_double,  &
            dim1name='pft', long_name='max height attained by a crop during year', units='m', &
            interpinic_flag='interp', readvar=readvar, data=this%htmx_patch)

       call restartvar(ncid=ncid, flag=flag,  varname='peaklai', xtype=ncd_int,  &
            dim1name='pft', long_name='Flag if at max allowed LAI or not', &
            flag_values=(/0,1/), nvalid_range=(/0,1/), &
            flag_meanings=(/'NOT-at-peak', 'AT_peak-LAI' /) , &
            interpinic_flag='interp', readvar=readvar, data=this%peaklai_patch)

       call restartvar(ncid=ncid, flag=flag,  varname='idop', xtype=ncd_int,  &
            dim1name='pft', long_name='Date of planting', units='jday', nvalid_range=(/1,366/), &
            interpinic_flag='interp', readvar=readvar, data=this%idop_patch)

       call restartvar(ncid=ncid, flag=flag,  varname='aleaf', xtype=ncd_double,  &
            dim1name='pft', long_name='leaf allocation coefficient', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%aleaf_patch)

       call restartvar(ncid=ncid, flag=flag,  varname='aleafi', xtype=ncd_double,  &
            dim1name='pft', long_name='Saved leaf allocation coefficient from phase 2', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%aleafi_patch)

       call restartvar(ncid=ncid, flag=flag,  varname='astem', xtype=ncd_double,  &
            dim1name='pft', long_name='stem allocation coefficient', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%astem_patch)

       call restartvar(ncid=ncid, flag=flag,  varname='astemi', xtype=ncd_double,  &
            dim1name='pft', long_name='Saved stem allocation coefficient from phase 2', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%astemi_patch)

       call restartvar(ncid=ncid, flag=flag,  varname='hdidx', xtype=ncd_double,  &
            dim1name='pft', long_name='cold hardening index', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%hdidx_patch)

       call restartvar(ncid=ncid, flag=flag,  varname='cumvd', xtype=ncd_double,  &
            dim1name='pft', long_name='cumulative vernalization d', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%cumvd_patch)

      call restartvar(ncid=ncid, flag=flag,  varname='gddmaturity', xtype=ncd_double,  &
            dim1name='pft', long_name='Growing degree days needed to harvest', units='ddays', &
            interpinic_flag='interp', readvar=readvar, data=this%gddmaturity_patch)

       call restartvar(ncid=ncid, flag=flag,  varname='huileaf', xtype=ncd_double,  &
            dim1name='pft', long_name='heat unit index needed from planting to leaf emergence', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%huileaf_patch)

       call restartvar(ncid=ncid, flag=flag,  varname='huigrain', xtype=ncd_double,  &
            dim1name='pft', long_name='heat unit index needed to reach vegetative maturity', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%huigrain_patch)

       call restartvar(ncid=ncid, flag=flag, varname='grain_flag', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%grain_flag_patch)

       ! Read or write variable(s) with mxharvests dimension
       ! BACKWARDS_COMPATIBILITY(ssr, 2022-03-31) See note in CallRestartvarDimOK()
       if (CallRestartvarDimOK(ncid, flag, 'mxharvests')) then
          call restartvar(ncid=ncid, flag=flag, varname='gddmaturity_thisyr', xtype=ncd_double,  &
               dim1name='pft', dim2name='mxharvests', switchdim=.true., &
               long_name='crop harvest dates for this patch this year', units='day of year', &
               scale_by_thickness=.false., &
               interpinic_flag='interp', readvar=readvar, data=this%gddmaturity_thisyr)
       end if
    end if
    if ( flag == 'read' .and. num_reseed_patch > 0 )then
       if ( masterproc ) write(iulog, *) 'Reseed dead plants for CNVegState'
       do i = 1, num_reseed_patch
          p = filter_reseed_patch(i)
          ! phenology variables
          this%dormant_flag_patch(p)   = 1._r8
          this%days_active_patch(p)    = 0._r8
          this%onset_flag_patch(p)     = 0._r8
          this%onset_counter_patch(p)  = 0._r8
          this%onset_gddflag_patch(p)  = 0._r8
          this%onset_fdd_patch(p)      = 0._r8
          this%onset_gdd_patch(p)      = 0._r8
          this%onset_swi_patch(p)      = 0._r8
          this%offset_flag_patch(p)    = 0._r8
          this%offset_counter_patch(p) = 0._r8
          this%offset_fdd_patch(p)     = 0._r8
          this%offset_swi_patch(p)     = 0._r8
          this%lgsf_patch(p)           = 0._r8
          this%bglfr_patch(p)          = 0._r8
          this%bgtr_patch(p)           = 0._r8
          this%annavg_t2m_patch(p)     = 280._r8
          this%tempavg_t2m_patch(p)    = 0._r8
          this%grain_flag_patch(p)     = 0._r8

          this%c_allometry_patch(p)           = 0._r8
          this%n_allometry_patch(p)           = 0._r8
          this%tempsum_potential_gpp_patch(p) = 0._r8
          this%annsum_potential_gpp_patch(p)  = 0._r8
          this%tempmax_retransn_patch(p)      = 0._r8
          this%annmax_retransn_patch(p)       = 0._r8
          this%downreg_patch(p)               = 0._r8
          this%leafcn_offset_patch(p) = spval
          this%plantCN_patch(p)       = spval
       end do
    end if

  end subroutine Restart

end module CNVegStateType
