module OzoneMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Calculates ozone-induced stress.
  !
  ! Note that the ozone calculations need to happen AFTER rssun and rsshade are computed
  ! by the Photosynthesis routine. However, Photosynthesis also uses the ozone stress
  ! computed here. Thus, the ozone stress computed in timestep i is applied in timestep
  ! (i+1), requiring these stresses to be saved on the restart file.
  !
  ! Developed by Danica Lombardozzi: Lombardozzi, D., S. Levis, G. Bonan, P. G. Hess, and
  ! J. P. Sparks (2015), The Influence of Chronic Ozone Exposure on Global Carbon and
  ! Water Cycles, J Climate, 28(1), 292–305, doi:10.1175/JCLI-D-14-00223.1.
  !
  ! !USES:
#include "shr_assert.h"
  use shr_kind_mod, only : r8 => shr_kind_r8
  use decompMod   , only : bounds_type
  use clm_varcon  , only : spval
  use clm_varctl  , only : iulog
  use OzoneBaseMod, only : ozone_base_type
  use abortutils  , only : endrun
  use PatchType   , only : patch
  use pftconMod   , only : pftcon

  implicit none
  save
  private

  ! !PUBLIC TYPES:
  type, extends(ozone_base_type), public :: ozone_type
     private
     ! Private data members
     integer :: stress_method  ! Which ozone stress parameterization we're using in this run

     real(r8), pointer :: o3uptakesha_patch(:) ! ozone dose, shaded leaves (mmol O3/m^2)
     real(r8), pointer :: o3uptakesun_patch(:) ! ozone dose, sunlit leaves (mmol O3/m^2)

     ! NOTE(wjs, 2014-09-29) tlai_old_patch really belongs alongside tlai_patch in
     ! CanopyStateType.  But there are problems with any way I can think to implement
     ! that:
     !
     ! - Updating tlai_old from a call in clm_driver, just before tlai is updated: This
     !   is problematic to do correctly because tlai is updated in different places
     !   depending on whether you're using SP, CN or ED.
     !
     ! - Updating tlai_old within each routine that updates tlai: This feels fragile,
     !   since it depends on each scheme remembering to do this update at the correct
     !   time.
     !
     ! - Making tlai a private member of CanopyFluxes, with getter and setter methods.
     !   Then the setter method would also set tlai_old. This feels like the most robust
     !   solution, but we don't have any precedent for using getters and setters for data
     !   arrays.
     real(r8), pointer :: tlai_old_patch(:)  ! tlai from last time step

   contains
     ! Public routines
     procedure, public :: Init
     procedure, public :: Restart
     procedure, public :: CalcOzoneUptake
     procedure, public :: CalcOzoneStress

     ! Private routines
     procedure, private :: InitAllocate
     procedure, private :: InitHistory
     procedure, private :: InitCold

     ! Calculate ozone uptake for a single point, for just sunlit or shaded leaves
     procedure, private, nopass :: CalcOzoneUptakeOnePoint

     ! Original ozone stress functions from Danica Lombardozzi 2015
     procedure, private         :: CalcOzoneStressLombardozzi2015           ! Stress parameterization
     procedure, private, nopass :: CalcOzoneStressLombardozzi2015OnePoint   ! Ozone stress calculation for single point

     ! Ozone stress functions from Stefanie Falk
     procedure, private         :: CalcOzoneStressFalk          ! Stress parameterization
     procedure, private, nopass :: CalcOzoneStressFalkOnePoint  ! Ozone stress calculation for single point

  end type ozone_type

  ! !PRIVATE TYPES:
  integer, parameter :: stress_method_lombardozzi2015 = 1
  integer, parameter :: stress_method_falk = 2

  ! TODO(wjs, 2014-09-29) The following parameters should eventually be moved to the
  ! params file. Parameters differentiated on veg type should be put on the params file
  ! with a pft dimension.

  ! o3:h2o resistance ratio defined by Sitch et al. 2007
  real(r8), parameter :: ko3 = 1.67_r8

  ! LAI threshold for LAIs that asymptote and don't reach 0
  real(r8), parameter :: lai_thresh = 0.5_r8

  ! threshold below which o3flux is set to 0 (nmol m^-2 s^-1)
  real(r8), parameter :: o3_flux_threshold = 0.8_r8

  ! o3 intercepts and slopes for photosynthesis
  real(r8), parameter :: needleleafPhotoInt   = 0.8390_r8  ! units = unitless
  real(r8), parameter :: needleleafPhotoSlope = 0._r8      ! units = per mmol m^-2
  real(r8), parameter :: broadleafPhotoInt    = 0.8752_r8  ! units = unitless
  real(r8), parameter :: broadleafPhotoSlope  = 0._r8      ! units = per mmol m^-2
  real(r8), parameter :: nonwoodyPhotoInt     = 0.8021_r8  ! units = unitless
  real(r8), parameter :: nonwoodyPhotoSlope   = -0.0009_r8 ! units = per mmol m^-2

  ! o3 intercepts and slopes for conductance
  real(r8), parameter :: needleleafCondInt    = 0.7823_r8  ! units = unitless
  real(r8), parameter :: needleleafCondSlope  = 0.0048_r8  ! units = per mmol m^-2
  real(r8), parameter :: broadleafCondInt     = 0.9125_r8  ! units = unitless
  real(r8), parameter :: broadleafCondSlope   = 0._r8      ! units = per mmol m^-2
  real(r8), parameter :: nonwoodyCondInt      = 0.7511_r8  ! units = unitless
  real(r8), parameter :: nonwoodyCondSlope    = 0._r8      ! units = per mmol m^-2

  ! Data is currently only available for broadleaf species (Dec 2020)
  ! o3 intercepts and slopes for JmaxO3/Jmax0
  real(r8), parameter :: needleleafJmaxInt   = 1._r8                ! units = unitless
  real(r8), parameter :: needleleafJmaxSlope = 0._r8                ! units = per mmol m^-2
  real(r8), parameter :: broadleafJmaxInt    = 1._r8                ! units = unitless
  real(r8), parameter :: broadleafJmaxSlope  = -0.0037_r8           ! units = per mmol m^-2
  real(r8), parameter :: nonwoodyJmaxInt     = 1._r8                ! units = unitless
  real(r8), parameter :: nonwoodyJmaxSlope   = 0._r8                ! units = per mmol m^-2


  character(len=*), parameter, private :: sourcefile = &
       __FILE__

contains

  ! ========================================================================
  ! Infrastructure routines (initialization, restart, etc.)
  ! ========================================================================

  !-----------------------------------------------------------------------
  subroutine Init(this, bounds, o3_veg_stress_method)
    !
    ! !DESCRIPTION:
    ! Initialize ozone data structure
    !
    ! !USES:
   use clm_varctl   , only : use_luna
    !
    ! !ARGUMENTS:
    class(ozone_type), intent(inout) :: this
    type(bounds_type), intent(in)    :: bounds
    character(len=*),  intent(in)    :: o3_veg_stress_method
    !-----------------------------------------------------------------------

    if (o3_veg_stress_method=='stress_lombardozzi2015') then
       this%stress_method = stress_method_lombardozzi2015
    else if (o3_veg_stress_method=='stress_falk') then
       this%stress_method = stress_method_falk
       if (.not. use_luna ) call endrun(' use_luna=.true. is required when o3_veg_stress_method = stress_falk.')
    else
       call endrun('unknown ozone stress method')
    end if

    call this%InitAllocate(bounds)
    call this%InitHistory(bounds)
    call this%InitCold(bounds)

  end subroutine Init


  !-----------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! !DESCRIPTION:
    ! Allocate memory for ozone data structure
    !
    ! !USES:
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
    !
    ! !ARGUMENTS:
    class(ozone_type), intent(inout) :: this
    type(bounds_type), intent(in)    :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    !-----------------------------------------------------------------------

    begp = bounds%begp
    endp = bounds%endp

    call this%InitAllocateBase(bounds)

    allocate(this%o3uptakesha_patch(begp:endp)) ; this%o3uptakesha_patch(:) = nan
    allocate(this%o3uptakesun_patch(begp:endp)) ; this%o3uptakesun_patch(:) = nan
    allocate(this%tlai_old_patch(begp:endp))    ; this%tlai_old_patch(:) = nan

  end subroutine InitAllocate

  !-----------------------------------------------------------------------
  subroutine InitHistory(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize ozone history variables
    !
    ! !USES:
    use histFileMod  , only : hist_addfld1d
    !
    ! !ARGUMENTS:
    class(ozone_type), intent(inout) :: this
    type(bounds_type), intent(in)    :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp

    character(len=*), parameter :: subname = 'InitHistory'
    !-----------------------------------------------------------------------

    begp = bounds%begp
    endp = bounds%endp

    this%o3uptakesun_patch(begp:endp) = spval
    call hist_addfld1d (fname='O3UPTAKESUN', units='mmol/m^2', &
         avgflag='A', long_name='total ozone flux into sunlit leaves', &
         ptr_patch=this%o3uptakesun_patch)

    this%o3uptakesha_patch(begp:endp) = spval
    call hist_addfld1d (fname='O3UPTAKESHA', units='mmol/m^2', &
         avgflag='A', long_name='total ozone flux into shaded leaves', &
         ptr_patch=this%o3uptakesha_patch)

    if (this%stress_method==stress_method_lombardozzi2015) then
       ! For this and the following variables: how should we include leaf area in the
       ! patch averaging?
       !
       ! Note that for now we need to reset these o3 coefficients over non-exposed veg
       ! points each time step for the sake of these diagnostics - to avoid weirdness and
       ! non-exact restarts due to coefficients over non-exposed patches "remembering"
       ! their value from when they were last exposed. (This resetting to 1 is done below,
       ! in the CalcOzoneStress* routines.) But if we could set the scaling so that the
       ! coefficients are only averaged over currently-exposed veg patches, then we could
       ! avoid needing to reset the coefficients to 1 each time step over non-exposed
       ! patches.
       !
       ! An alternative would be to set these values to spval over non-exposed patches:
       ! that way, averages would just be taken over exposed patches (as opposed to
       ! averaging in 1 values for non-exposed patches). However, that would conflict with
       ! the goals of https://github.com/ESCOMP/CTSM/projects/35 (specifically, having
       ! frac fields that can be used to determine the appropriate weighting of each grid
       ! cell in producing regional / global averages), so if we want to exclude
       ! non-exposed patches from the gridcell average, we should probably come up with a
       ! way to do so via a more explicit mechanism similar to l2g_scale_type (maybe with
       ! a new p2c_scale_type that allows excluding non-exposed patches???).
       this%o3coefvsha_patch(begp:endp) = spval
       call hist_addfld1d (fname='O3COEFPSNSHA', units='unitless', &
            avgflag='A', long_name='ozone coefficient for photosynthesis for shaded leaves', &
            ptr_patch=this%o3coefvsha_patch, l2g_scale_type='veg')

       this%o3coefvsun_patch(begp:endp) = spval
       call hist_addfld1d (fname='O3COEFPSNSUN', units='unitless', &
            avgflag='A', long_name='ozone coefficient for photosynthesis for sunlit leaves', &
            ptr_patch=this%o3coefvsun_patch, l2g_scale_type='veg')

       this%o3coefgsha_patch(begp:endp) = spval
       call hist_addfld1d (fname='O3COEFGSSHA', units='unitless', &
            avgflag='A', long_name='ozone coefficient for stomatal conductance for shaded leaves', &
            ptr_patch=this%o3coefgsha_patch, l2g_scale_type='veg')

       this%o3coefgsun_patch(begp:endp) = spval
       call hist_addfld1d (fname='O3COEFGSSUN', units='unitless', &
            avgflag='A', long_name='ozone coefficient for stomatal conductance for sunlit leaves', &
            ptr_patch=this%o3coefgsun_patch, l2g_scale_type='veg')

    elseif (this%stress_method==stress_method_falk) then
       !
       this%o3coefjmaxsha_patch(begp:endp) = spval
       call hist_addfld1d (fname='O3COEFJMAXSHA', units='unitless', &
            avgflag='A', long_name='ozone coefficient for maximum electron transport rate for shaded leaves', &
            ptr_patch=this%o3coefjmaxsha_patch, l2g_scale_type='veg')

       this%o3coefjmaxsun_patch(begp:endp) = spval
       call hist_addfld1d (fname='O3COEFJMAXSUN', units='unitless', &
            avgflag='A', long_name='ozone coefficient for maximum electron transport rate sunlit leaves', &
            ptr_patch=this%o3coefjmaxsun_patch, l2g_scale_type='veg')
    else
       call endrun('unknown ozone stress method')
    end if

  end subroutine InitHistory

  !-----------------------------------------------------------------------
  subroutine InitCold(this, bounds)
    !
    ! !DESCRIPTION:
    ! Perform cold-start initialization for ozone
    !
    ! !ARGUMENTS:
    class(ozone_type), intent(inout) :: this
    type(bounds_type), intent(in)    :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp

    character(len=*), parameter :: subname = 'InitCold'
    !-----------------------------------------------------------------------

    begp = bounds%begp
    endp = bounds%endp

    call this%InitColdBase(bounds)

    this%o3uptakesha_patch(begp:endp) = 0._r8
    this%o3uptakesun_patch(begp:endp) = 0._r8
    this%tlai_old_patch(begp:endp) = 0._r8

  end subroutine InitCold

  !-----------------------------------------------------------------------
  subroutine Restart(this, bounds, ncid, flag)
    !
    ! !DESCRIPTION:
    ! Handle restart of ozone variables.
    !
    ! !USES:
    use ncdio_pio  , only : file_desc_t, ncd_inqvdlen, ncd_double
    use restUtilMod
    !
    ! !ARGUMENTS:
    class(ozone_type) :: this
    type(bounds_type), intent(in)    :: bounds
    type(file_desc_t), intent(inout) :: ncid   ! netcdf id
    character(len=*) , intent(in)    :: flag   ! 'read', 'write' or 'define'
    !
    ! !LOCAL VARIABLES:
    logical :: readvar

    character(len=*), parameter :: subname = 'Restart'
    !-----------------------------------------------------------------------

    call restartvar(ncid=ncid, flag=flag, varname='o3_tlaiold', xtype=ncd_double, &
         dim1name='pft', &
         long_name='one-sided leaf area index, from previous timestep, for ozone calculations', units='', &
         readvar=readvar, interpinic_flag='interp', data=this%tlai_old_patch)

    call restartvar(ncid=ncid, flag=flag, varname='o3uptakesha', xtype=ncd_double, &
         dim1name='pft', &
         long_name='ozone uptake for shaded leaves', units='mmol m^-3', &
         readvar=readvar, interpinic_flag='interp', data=this%o3uptakesha_patch)

    call restartvar(ncid=ncid, flag=flag, varname='o3uptakesun', xtype=ncd_double, &
         dim1name='pft', &
         long_name='ozone uptake for sunlit leaves', units='mmol m^-3', &
         readvar=readvar, interpinic_flag='interp', data=this%o3uptakesun_patch)

  end subroutine Restart

  ! ========================================================================
  ! Science routines
  ! ========================================================================

  !-----------------------------------------------------------------------
  subroutine CalcOzoneUptake(this, bounds, num_exposedvegp, filter_exposedvegp, &
       forc_pbot, forc_th, rssun, rssha, rb, ram, tlai, forc_o3)
    !
    ! !DESCRIPTION:
    ! Calculate ozone uptake.
    !
    ! !ARGUMENTS:
    class(ozone_type)      , intent(inout) :: this
    type(bounds_type)      , intent(in)    :: bounds
    integer  , intent(in) :: num_exposedvegp           ! number of points in filter_exposedvegp
    integer  , intent(in) :: filter_exposedvegp(:)     ! patch filter for non-snow-covered veg
    real(r8) , intent(in) :: forc_pbot( bounds%begc: ) ! atmospheric pressure (Pa)
    real(r8) , intent(in) :: forc_th( bounds%begc: )   ! atmospheric potential temperature (K)
    real(r8) , intent(in) :: rssun( bounds%begp: )     ! leaf stomatal resistance, sunlit leaves (s/m)
    real(r8) , intent(in) :: rssha( bounds%begp: )     ! leaf stomatal resistance, shaded leaves (s/m)
    real(r8) , intent(in) :: rb( bounds%begp: )        ! boundary layer resistance (s/m)
    real(r8) , intent(in) :: ram( bounds%begp: )       ! aerodynamical resistance (s/m)
    real(r8) , intent(in) :: tlai( bounds%begp: )      ! one-sided leaf area index, no burying by snow
    real(r8) , intent(in) :: forc_o3( bounds%begg: )   ! ozone partial pressure (mol/mol)
    !
    ! !LOCAL VARIABLES:
    integer  :: fp             ! filter index
    integer  :: p              ! patch index
    integer  :: c              ! column index
    integer  :: g              ! gridcell index

    character(len=*), parameter :: subname = 'CalcOzoneUptake'
    !-----------------------------------------------------------------------

    ! Enforce expected array sizes
    SHR_ASSERT_ALL_FL((ubound(forc_pbot) == (/bounds%endc/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(forc_th) == (/bounds%endc/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(forc_o3) == (/bounds%endg/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(rssun) == (/bounds%endp/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(rssha) == (/bounds%endp/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(rb) == (/bounds%endp/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(ram) == (/bounds%endp/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(tlai) == (/bounds%endp/)), sourcefile, __LINE__)

    associate( &
         o3uptakesha => this%o3uptakesha_patch                , & ! Output: [real(r8) (:)] ozone dose
         o3uptakesun => this%o3uptakesun_patch                , & ! Output: [real(r8) (:)] ozone dose
         tlai_old    => this%tlai_old_patch                     & ! Output: [real(r8) (:)] tlai from last time step
         )

      do fp = 1, num_exposedvegp
         p = filter_exposedvegp(fp)
         c = patch%column(p)
         g = patch%gridcell(p)

         ! Ozone uptake for shaded leaves
         call CalcOzoneUptakeOnePoint( &
              forc_ozone=forc_o3(g), forc_pbot=forc_pbot(c), forc_th=forc_th(c), &
              rs=rssha(p), rb=rb(p), ram=ram(p), &
              tlai=tlai(p), tlai_old=tlai_old(p), pft_type=patch%itype(p), &
              o3uptake=o3uptakesha(p))

         ! Ozone uptake for sunlit leaves
         call CalcOzoneUptakeOnePoint( &
              forc_ozone=forc_o3(g), forc_pbot=forc_pbot(c), forc_th=forc_th(c), &
              rs=rssun(p), rb=rb(p), ram=ram(p), &
              tlai=tlai(p), tlai_old=tlai_old(p), pft_type=patch%itype(p), &
              o3uptake=o3uptakesun(p))

         tlai_old(p) = tlai(p)

      end do

    end associate

  end subroutine CalcOzoneUptake

  !-----------------------------------------------------------------------
  subroutine CalcOzoneUptakeOnePoint( &
       forc_ozone, forc_pbot, forc_th, &
       rs, rb, ram, &
       tlai, tlai_old, pft_type, &
       o3uptake)
    !
    ! !DESCRIPTION:
    ! Calculates ozone uptake for a single point, for just sunlit or shaded leaves
    !
    ! !USES:
    use shr_const_mod        , only : SHR_CONST_RGAS
    use clm_time_manager     , only : get_step_size
    !
    ! !ARGUMENTS:
    real(r8) , intent(in)    :: forc_ozone ! ozone partial pressure (mol/mol)
    real(r8) , intent(in)    :: forc_pbot  ! atmospheric pressure (Pa)
    real(r8) , intent(in)    :: forc_th    ! atmospheric potential temperature (K)
    real(r8) , intent(in)    :: rs         ! leaf stomatal resistance (s/m)
    real(r8) , intent(in)    :: rb         ! boundary layer resistance (s/m)
    real(r8) , intent(in)    :: ram        ! aerodynamical resistance (s/m)
    real(r8) , intent(in)    :: tlai       ! one-sided leaf area index, no burying by snow
    real(r8) , intent(in)    :: tlai_old   ! tlai from last time step
    integer  , intent(in)    :: pft_type   ! vegetation type, for indexing into pftvarcon arrays
    real(r8) , intent(inout) :: o3uptake   ! ozone entering the leaf
    !
    ! !LOCAL VARIABLES:
    integer  :: dtime          ! land model time step (sec)
    real(r8) :: dtimeh         ! time step in hours
    real(r8) :: o3concnmolm3   ! o3 concentration (nmol/m^3)
    real(r8) :: o3flux         ! instantaneous o3 flux (nmol m^-2 s^-1)
    real(r8) :: o3fluxcrit     ! instantaneous o3 flux beyond threshold (nmol m^-2 s^-1)
    real(r8) :: o3fluxperdt    ! o3 flux per timestep (mmol m^-2)
    real(r8) :: heal           ! o3uptake healing rate based on % of new leaves growing (mmol m^-2)
    real(r8) :: leafturn       ! leaf turnover time / mortality rate (per hour)
    real(r8) :: decay          ! o3uptake decay rate based on leaf lifetime (mmol m^-2)

    character(len=*), parameter :: subname = 'CalcOzoneUptakeOnePoint'
    !-----------------------------------------------------------------------

    ! convert o3 from mol/mol to nmol m^-3
    o3concnmolm3 = forc_ozone * 1.e9_r8 * (forc_pbot/(forc_th*SHR_CONST_RGAS*0.001_r8))

    ! calculate instantaneous flux
    o3flux = o3concnmolm3/ (ko3*rs+ rb + ram)

    ! apply o3 flux threshold
    if (o3flux < o3_flux_threshold) then
       o3fluxcrit = 0._r8
    else
       o3fluxcrit = o3flux - o3_flux_threshold
    endif

    dtime  = get_step_size()
    dtimeh = dtime / 3600._r8

    ! calculate o3 flux per timestep
    o3fluxperdt = o3fluxcrit * dtime * 0.000001_r8

    if (tlai > lai_thresh) then
       ! checking if new leaf area was added
       if (tlai - tlai_old > 0) then
          ! minimizing o3 damage to new leaves
          heal = max(0._r8,(((tlai-tlai_old)/tlai)*o3fluxperdt))
       else
          heal = 0._r8
       endif

       if (pftcon%evergreen(pft_type) == 1) then
          leafturn = 1._r8/(pftcon%leaf_long(pft_type)*365._r8*24._r8)
       else
          leafturn = 0._r8
       endif

       ! o3 uptake decay based on leaf lifetime for evergreen plants
       decay = o3uptake * leafturn * dtimeh
       !cumulative uptake (mmol m^-2)
       o3uptake = max(0._r8, o3uptake + o3fluxperdt - decay - heal)

    else
       o3uptake = 0._r8
    end if

  end subroutine CalcOzoneUptakeOnePoint

  !-----------------------------------------------------------------------
  subroutine CalcOzoneStress(this, bounds, &
       num_exposedvegp, filter_exposedvegp, &
       num_noexposedvegp, filter_noexposedvegp)
    !
    ! !DESCRIPTION:
    ! Calculate ozone stress.
    !
    ! !ARGUMENTS:
    class(ozone_type), intent(inout) :: this
    type(bounds_type), intent(in) :: bounds
    integer , intent(in) :: num_exposedvegp         ! number of points in filter_exposedvegp
    integer , intent(in) :: filter_exposedvegp(:)   ! patch filter for non-snow-covered veg
    integer , intent(in) :: num_noexposedvegp       ! number of points in filter_noexposedvegp
    integer , intent(in) :: filter_noexposedvegp(:) ! patch filter for veg where frac_veg_nosno is 0
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'CalcOzoneStress'
    !-----------------------------------------------------------------------

    select case (this%stress_method)
    case (stress_method_lombardozzi2015)
       call this%CalcOzoneStressLombardozzi2015(bounds, &
            num_exposedvegp, filter_exposedvegp, &
            num_noexposedvegp, filter_noexposedvegp)
    case (stress_method_falk)
       call this%CalcOzoneStressFalk(bounds, &
            num_exposedvegp, filter_exposedvegp, &
            num_noexposedvegp, filter_noexposedvegp)
    case default
       write(iulog,*) 'ERROR: unknown ozone stress method: ', this%stress_method
       call endrun('Unknown ozone stress method')
    end select

  end subroutine CalcOzoneStress

  !-----------------------------------------------------------------------
  subroutine CalcOzoneStressLombardozzi2015(this, bounds, &
       num_exposedvegp, filter_exposedvegp, &
       num_noexposedvegp, filter_noexposedvegp)
    !
    ! !DESCRIPTION:
    ! Calculate ozone stress.
    !
    ! This subroutine uses the Lombardozzi2015 formulation for ozone stress
    !
    ! !ARGUMENTS:
    class(ozone_type), intent(inout) :: this
    type(bounds_type), intent(in) :: bounds
    integer , intent(in) :: num_exposedvegp         ! number of points in filter_exposedvegp
    integer , intent(in) :: filter_exposedvegp(:)   ! patch filter for non-snow-covered veg
    integer , intent(in) :: num_noexposedvegp       ! number of points in filter_noexposedvegp
    integer , intent(in) :: filter_noexposedvegp(:) ! patch filter for veg where frac_veg_nosno is 0
    !
    ! !LOCAL VARIABLES:
    integer  :: fp             ! filter index
    integer  :: p              ! patch index

    character(len=*), parameter :: subname = 'CalcOzoneStressLombardozzi2015'
    !-----------------------------------------------------------------------

    associate( &
         o3uptakesha => this%o3uptakesha_patch                , & ! Input:  [real(r8) (:)] ozone dose
         o3uptakesun => this%o3uptakesun_patch                , & ! Input:  [real(r8) (:)] ozone dose
         o3coefvsha  => this%o3coefvsha_patch                 , & ! Output: [real(r8) (:)] ozone coef
         o3coefvsun  => this%o3coefvsun_patch                 , & ! Output: [real(r8) (:)] ozone coef
         o3coefgsha  => this%o3coefgsha_patch                 , & ! Output: [real(r8) (:)] ozone coef
         o3coefgsun  => this%o3coefgsun_patch                   & ! Output: [real(r8) (:)] ozone coef
         )

      do fp = 1, num_exposedvegp
         p = filter_exposedvegp(fp)

         ! Ozone stress for shaded leaves
         call CalcOzoneStressLombardozzi2015OnePoint( &
              pft_type=patch%itype(p), o3uptake=o3uptakesha(p), &
              o3coefv=o3coefvsha(p), o3coefg=o3coefgsha(p))

         ! Ozone stress for sunlit leaves
         call CalcOzoneStressLombardozzi2015OnePoint( &
              pft_type=patch%itype(p), o3uptake=o3uptakesun(p), &
              o3coefv=o3coefvsun(p), o3coefg=o3coefgsun(p))
      end do

      do fp = 1, num_noexposedvegp
         p = filter_noexposedvegp(fp)

         ! See notes above in InitHistory about why these need to be set to 1 over
         ! non-exposed veg points each time step.
         o3coefvsha(p) = 1._r8
         o3coefgsha(p) = 1._r8
         o3coefvsun(p) = 1._r8
         o3coefgsun(p) = 1._r8
      end do

    end associate

  end subroutine CalcOzoneStressLombardozzi2015

  !-----------------------------------------------------------------------
  subroutine CalcOzoneStressLombardozzi2015OnePoint( &
       pft_type, o3uptake, &
       o3coefv, o3coefg)
    !
    ! !DESCRIPTION:
    ! Calculates ozone stress for a single point, for just sunlit or shaded leaves
    !
    ! This subroutine uses the Lombardozzi2015 formulation for ozone stress
    !
    ! !ARGUMENTS:
    integer  , intent(in)    :: pft_type   ! vegetation type, for indexing into pftvarcon arrays
    real(r8) , intent(in)    :: o3uptake   ! ozone entering the leaf
    real(r8) , intent(out)   :: o3coefv    ! ozone coefficient for photosynthesis (0 - 1)
    real(r8) , intent(out)   :: o3coefg    ! ozone coefficient for conductance (0 - 1)
    !
    ! !LOCAL VARIABLES:
    real(r8) :: photoInt       ! intercept for photosynthesis
    real(r8) :: photoSlope     ! slope for photosynthesis
    real(r8) :: condInt        ! intercept for conductance
    real(r8) :: condSlope      ! slope for conductance

    character(len=*), parameter :: subname = 'CalcOzoneStressLombardozzi2015OnePoint'
    !-----------------------------------------------------------------------

    if (o3uptake == 0._r8) then
       ! No o3 damage if no o3 uptake
       o3coefv = 1._r8
       o3coefg = 1._r8
    else
       ! Determine parameter values for this pft
       ! TODO(wjs, 2014-10-01) Once these parameters are moved into the params file, this
       ! logic can be removed.
       if (pft_type>3) then
          if (pftcon%woody(pft_type)==0) then
             photoInt   = nonwoodyPhotoInt
             photoSlope = nonwoodyPhotoSlope
             condInt    = nonwoodyCondInt
             condSlope  = nonwoodyCondSlope
          else
             photoInt   = broadleafPhotoInt
             photoSlope = broadleafPhotoSlope
             condInt    = broadleafCondInt
             condSlope  = broadleafCondSlope
          end if
       else
          photoInt   = needleleafPhotoInt
          photoSlope = needleleafPhotoSlope
          condInt    = needleleafCondInt
          condSlope  = needleleafCondSlope
       end if

       ! Apply parameter values to compute o3 coefficients
       o3coefv = max(0._r8, min(1._r8, photoInt + photoSlope * o3uptake))
       o3coefg = max(0._r8, min(1._r8, condInt  + condSlope  * o3uptake))

    end if

  end subroutine CalcOzoneStressLombardozzi2015OnePoint


  !-----------------------------------------------------------------------
  subroutine CalcOzoneStressFalk(this, bounds, &
       num_exposedvegp, filter_exposedvegp, &
       num_noexposedvegp, filter_noexposedvegp)
    !
    ! !DESCRIPTION:
    ! Calculate ozone stress.
    !
    ! This subroutine uses the Falk formulation for ozone stress
    !
    ! !USES:
    use LunaMod        , only : is_time_to_run_LUNA
    !
    ! !ARGUMENTS:
    class(ozone_type), intent(inout) :: this
    type(bounds_type), intent(in)    :: bounds
    integer , intent(in) :: num_exposedvegp         ! number of points in filter_exposedvegp
    integer , intent(in) :: filter_exposedvegp(:)   ! patch filter for non-snow-covered veg
    integer , intent(in) :: num_noexposedvegp       ! number of points in filter_noexposedvegp
    integer , intent(in) :: filter_noexposedvegp(:) ! patch filter for veg where frac_veg_nosno is 0
    !
    ! !LOCAL VARIABLES:
    integer  :: fp             ! filter index
    integer  :: p              ! patch index

    character(len=*), parameter :: subname = 'CalcOzoneStressFalk'
    !-----------------------------------------------------------------------

    if (is_time_to_run_LUNA()) then

      associate( &
         o3uptakesha => this%o3uptakesha_patch                 , & ! Input:  [real(r8) (:)] ozone dose
         o3uptakesun => this%o3uptakesun_patch                 , & ! Input:  [real(r8) (:)] ozone dose
         o3coefjmaxsha => this%o3coefjmaxsha_patch             , & ! Output: [real(r8) (:)] ozone coef jmax sha
         o3coefjmaxsun => this%o3coefjmaxsun_patch               & ! Output: [real(r8) (:)] ozone coef jmax sun
         )

      do fp = 1, num_exposedvegp
         p = filter_exposedvegp(fp)

         ! Ozone stress for shaded leaves
          call CalcOzoneStressFalkOnePoint( &
            pft_type=patch%itype(p), o3uptake=o3uptakesha(p), &
            o3coefjmax=o3coefjmaxsha(p))

         ! Ozone stress for sunlit leaves
         call CalcOzoneStressFalkOnePoint( &
            pft_type=patch%itype(p), o3uptake=o3uptakesun(p), &
            o3coefjmax=o3coefjmaxsun(p))

      end do

      do fp = 1, num_noexposedvegp
         p = filter_noexposedvegp(fp)

         ! See notes above in InitHistory about why these need to be set to 1 over
         ! non-exposed veg points each time step.
         o3coefjmaxsha(p) = 1._r8
         o3coefjmaxsun(p) = 1._r8
      end do

      end associate

   end if

  end subroutine CalcOzoneStressFalk

  !-----------------------------------------------------------------------
  subroutine CalcOzoneStressFalkOnePoint( pft_type, o3uptake, o3coefjmax)
    !
    ! !DESCRIPTION:
    ! Calculates ozone stress for a single point, for just sunlit or shaded leaves
    !
    ! This subroutine uses the Falk formulation for ozone stress
    !
    ! !ARGUMENTS:
    integer  , intent(in)        :: pft_type   ! vegetation type, for indexing into pftvarcon arrays
    real(r8) , intent(in)        :: o3uptake   ! ozone entering the leaf
    real(r8) , intent(inout)     :: o3coefjmax ! ozone coefficient for max. electron transport rate
    !
    ! !LOCAL VARIABLES:
    real(r8) :: jmaxInt        ! intercept for max. electron transport rate
    real(r8) :: jmaxSlope      ! slope for max. electron transport rate
    character(len=*), parameter :: subname = 'CalcOzoneStressFalkOnePoint'
    !-----------------------------------------------------------------------

    if (o3uptake == 0._r8) then
       o3coefjmax = 1._r8
    else
       ! Determine parameter values for this pft
       ! TODO(wjs, 2014-10-01) Once these parameters are moved into the params file, this
       ! logic can be removed.
       if (pft_type>3) then
          if (pftcon%woody(pft_type)==0) then
             jmaxInt    = nonwoodyJmaxInt
             jmaxSlope  = nonwoodyJmaxSlope
          else
             jmaxInt    = broadleafJmaxInt
             jmaxSlope  = broadleafJmaxSlope
          end if
       else
          jmaxInt    = needleleafJmaxInt
          jmaxSlope  = needleleafJmaxSlope
       end if
       ! Apply parameter values to compute o3 coefficients
       o3coefjmax = max(0._r8, min(1._r8, jmaxInt  + jmaxSlope  * o3uptake))
    end if

  end subroutine CalcOzoneStressFalkOnePoint

end module OzoneMod
