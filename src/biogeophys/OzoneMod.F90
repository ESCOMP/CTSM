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
  ! Developed by Danica Lombardozzi.
  !
  ! !USES:
#include "shr_assert.h"
  use shr_kind_mod, only : r8 => shr_kind_r8
  use decompMod   , only : bounds_type
  use clm_varcon  , only : spval
  use shr_log_mod , only : errMsg => shr_log_errMsg
  use OzoneBaseMod, only : ozone_base_type
  use abortutils  , only : endrun

  implicit none
  save
  private

  ! !PUBLIC TYPES:
  type, extends(ozone_base_type), public :: ozone_type
     private
     ! Private data members
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
     procedure, public :: CalcOzoneStress

     ! Private routines
     procedure, private :: InitAllocate
     procedure, private :: InitHistory
     procedure, private :: InitCold

     ! Calculate ozone stress for a single point, for just sunlit or shaded leaves
     procedure, private, nopass :: CalcOzoneStressOnePoint
  end type ozone_type

  interface ozone_type
     module procedure constructor
  end interface ozone_type

  ! !PRIVATE TYPES:
  
  ! TODO(wjs, 2014-09-29) This parameter will eventually become a spatially-varying
  ! value, obtained from ATM
  real(r8), parameter :: forc_ozone = 100._r8 * 1.e-9_r8  ! ozone partial pressure [mol/mol]

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

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

contains
  
  ! ========================================================================
  ! Infrastructure routines (initialization, restart, etc.)
  ! ========================================================================

  !-----------------------------------------------------------------------
  function constructor() result(ozone)
    !
    ! !DESCRIPTION:
    ! Return an instance of ozone_type
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(ozone_type) :: ozone  ! function result
    !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'constructor'
    !-----------------------------------------------------------------------

    ! DO NOTHING (simply return a variable of the appropriate type)

    ! Eventually this should call the Init routine (or replace the Init routine
    ! entirely). But I think it would be confusing to do that until we switch everything
    ! to use a constructor rather than the init routine.
    
  end function constructor


  !-----------------------------------------------------------------------
  subroutine Init(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize ozone data structure
    !
    ! !ARGUMENTS:
    class(ozone_type), intent(inout) :: this
    type(bounds_type), intent(in)    :: bounds
    !-----------------------------------------------------------------------

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

    call restartvar(ncid=ncid, flag=flag, varname='o3coefvsun', xtype=ncd_double, &
         dim1name='pft', &
         long_name='ozone coefficient for photosynthesis for sunlit leaves', units='unitless', &
         readvar=readvar, interpinic_flag='interp', data=this%o3coefvsun_patch)

    call restartvar(ncid=ncid, flag=flag, varname='o3coefgsun', xtype=ncd_double, &
         dim1name='pft', &
         long_name='ozone coefficient for stomatal conductance for sunlit leaves', units='unitless', &
         readvar=readvar, interpinic_flag='interp', data=this%o3coefgsun_patch)

    call restartvar(ncid=ncid, flag=flag, varname='o3coefvsha', xtype=ncd_double, &
         dim1name='pft', &
         long_name='ozone coefficient for photosynthesis for shaded leaves', units='unitless', &
         readvar=readvar, interpinic_flag='interp', data=this%o3coefvsha_patch)

    call restartvar(ncid=ncid, flag=flag, varname='o3coefgsha', xtype=ncd_double, &
         dim1name='pft', &
         long_name='ozone coefficient for stomatal conductance for shaded leaves', units='unitless', &
         readvar=readvar, interpinic_flag='interp', data=this%o3coefgsha_patch)

  end subroutine Restart

  ! ========================================================================
  ! Science routines
  ! ========================================================================

  !-----------------------------------------------------------------------
  subroutine CalcOzoneStress(this, bounds, num_exposedvegp, filter_exposedvegp, &
       forc_pbot, forc_th, rssun, rssha, rb, ram, tlai)
    !
    ! !DESCRIPTION:
    ! Calculate ozone stress.
    !
    ! !USES:
    use PatchType            , only : patch
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
    !
    ! !LOCAL VARIABLES:
    integer  :: fp             ! filter index
    integer  :: p              ! patch index
    integer  :: c              ! column index

    character(len=*), parameter :: subname = 'CalcOzoneStress'
    !-----------------------------------------------------------------------
    
    ! Enforce expected array sizes
    SHR_ASSERT_ALL((ubound(forc_pbot) == (/bounds%endc/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(forc_th) == (/bounds%endc/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(rssun) == (/bounds%endp/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(rssha) == (/bounds%endp/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(rb) == (/bounds%endp/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(ram) == (/bounds%endp/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(tlai) == (/bounds%endp/)), errMsg(sourcefile, __LINE__))

    associate( &
         o3coefvsha  => this%o3coefvsha_patch                 , & ! Output: [real(r8) (:)] ozone coef
         o3coefvsun  => this%o3coefvsun_patch                 , & ! Output: [real(r8) (:)] ozone coef
         o3coefgsha  => this%o3coefgsha_patch                 , & ! Output: [real(r8) (:)] ozone coef
         o3coefgsun  => this%o3coefgsun_patch                 , & ! Output: [real(r8) (:)] ozone coef
         o3uptakesha => this%o3uptakesha_patch                , & ! Output: [real(r8) (:)] ozone dose
         o3uptakesun => this%o3uptakesun_patch                , & ! Output: [real(r8) (:)] ozone dose
         tlai_old    => this%tlai_old_patch                     & ! Output: [real(r8) (:)] tlai from last time step
         )
         
    do fp = 1, num_exposedvegp
       p = filter_exposedvegp(fp)
       c = patch%column(p)

       ! Ozone stress for shaded leaves
       call CalcOzoneStressOnePoint( &
            forc_ozone=forc_ozone, forc_pbot=forc_pbot(c), forc_th=forc_th(c), &
            rs=rssha(p), rb=rb(p), ram=ram(p), &
            tlai=tlai(p), tlai_old=tlai_old(p), pft_type=patch%itype(p), &
            o3uptake=o3uptakesha(p), o3coefv=o3coefvsha(p), o3coefg=o3coefgsha(p))

       ! Ozone stress for sunlit leaves
       call CalcOzoneStressOnePoint( &
            forc_ozone=forc_ozone, forc_pbot=forc_pbot(c), forc_th=forc_th(c), &
            rs=rssun(p), rb=rb(p), ram=ram(p), &
            tlai=tlai(p), tlai_old=tlai_old(p), pft_type=patch%itype(p), &
            o3uptake=o3uptakesun(p), o3coefv=o3coefvsun(p), o3coefg=o3coefgsun(p))

       tlai_old(p) = tlai(p)

    end do

    end associate

  end subroutine CalcOzoneStress

  !-----------------------------------------------------------------------
  subroutine CalcOzoneStressOnePoint( &
       forc_ozone, forc_pbot, forc_th, &
       rs, rb, ram, &
       tlai, tlai_old, pft_type, &
       o3uptake, o3coefv, o3coefg)
    !
    ! !DESCRIPTION:
    ! Calculates ozone stress for a single point, for just sunlit or shaded leaves
    !
    ! !USES:
    use shr_const_mod        , only : SHR_CONST_RGAS
    use pftconMod            , only : pftcon
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
    real(r8) , intent(out)   :: o3coefv    ! ozone coefficient for photosynthesis (0 - 1)
    real(r8) , intent(out)   :: o3coefg    ! ozone coefficient for conductance (0 - 1)
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
    real(r8) :: photoInt       ! intercept for photosynthesis
    real(r8) :: photoSlope     ! slope for photosynthesis
    real(r8) :: condInt        ! intercept for conductance
    real(r8) :: condSlope      ! slope for conductance

    character(len=*), parameter :: subname = 'CalcOzoneStressOnePoint'
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

  end subroutine CalcOzoneStressOnePoint


end module OzoneMod
