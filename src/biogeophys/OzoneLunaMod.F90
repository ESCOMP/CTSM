module OzoneLunaMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Calculates ozone-induced stress for use with LUNA model.
  !
  ! A correction factor (1->0) in response to cumulative ozone uptake (CUO) 
  ! is calculated for Jmax. We assume a linear relationship between Jmax(CUO)
  ! normalized to Jmax(0). This replaces the direkt correction on stomatal 
  ! conductance and net photosynthesis in OzoneMod.
  !
  ! Developed by Stefanie Falk.
  !
  ! !USES:
#include "shr_assert.h"
  use shr_kind_mod, only : r8 => shr_kind_r8
  use decompMod   , only : bounds_type
  use clm_varcon  , only : spval
  use shr_log_mod , only : errMsg => shr_log_errMsg
  use OzoneMod, only : ozone_type
  use abortutils  , only : endrun

  implicit none
  save
  private

  ! !PUBLIC TYPES:
  type, extends(ozone_type), public :: ozone_luna_type
     private
     ! Private data members are inherited from base classes
     
  
   contains
     ! Public routines
     procedure, public :: Init
     procedure, public :: Restart
     procedure, public :: CalcOzoneStress
     
     ! Private routines
     procedure, private :: InitAllocateLuna
     procedure, private :: InitHistoryLuna
     procedure, private :: InitColdLuna

     ! Calculate ozone stress for sunlit or shaded leaves accumulated over one day
     procedure, public :: Acc24_OzoneStress_Luna
     procedure, private, nopass :: Acc24_OzoneStressOnePoint_Luna
  end type ozone_luna_type

  interface ozone_luna_type
     module procedure constructor
  end interface ozone_luna_type

  ! !PRIVATE TYPES:
  
  ! TODO(wjs, 2014-09-29) The following parameters should eventually be moved to the
  ! params file. Parameters differentiated on veg type should be put on the params file
  ! with a pft dimension.

  ! o3:h2o resistance ratio defined by Sitch et al. 2007
  ! This should actually be defined directly via diffusivities in air DH2O : DO3
  real(r8), parameter :: ko3 = 1.67_r8

  ! Data is currently only available for broadleaf species (Dec 2020)
  ! O3 intercepts and slopes for JmaxO3/Jmax0
  real(r8), parameter :: needleleafJmaxInt   = 1._r8           ! units = unitless 
  real(r8), parameter :: needleleafJmaxSlope = 0._r8           ! units = per mmol m^-2
  real(r8), parameter :: broadleafJmaxInt    = 1._r8           ! units = unitless
  real(r8), parameter :: broadleafJmaxSlope  = -0.0037_r8      ! units = per mmol m^-2
  real(r8), parameter :: nonwoodyJmaxInt     = 1._r8           ! units = unitless
  real(r8), parameter :: nonwoodyJmaxSlope   = 0._r8           ! units = per mmol m^-2
  
  character(len=*), parameter, private :: sourcefile = &
       __FILE__

contains
  ! ========================================================================
  ! Infrastructure routines (initialization, restart, etc.)
  ! ========================================================================

  !-----------------------------------------------------------------------
  function constructor() result(ozone_luna)
    !
    ! !DESCRIPTION:
    ! Return an instance of ozone_luna_type
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(ozone_luna_type) :: ozone_luna  ! function result
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
    class(ozone_luna_type), intent(inout) :: this
    type(bounds_type), intent(in)    :: bounds
    !-----------------------------------------------------------------------

    call this%InitAllocateLuna(bounds)
    call this%InitHistoryLuna(bounds)
    call this%InitColdLuna(bounds)
    
  end subroutine Init

!-----------------------------------------------------------------------
  subroutine InitAllocateLuna(this, bounds)
    !
    ! !DESCRIPTION:
    ! Allocate memory for ozone data structure
    !
    ! !USES:
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
    !
    ! !ARGUMENTS:
    class(ozone_luna_type), intent(inout) :: this
    type(bounds_type), intent(in)    :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    !-----------------------------------------------------------------------

    begp = bounds%begp
    endp = bounds%endp

    ! Allocating variables from OzoneMod
    call this%InitAllocate(bounds)

    ! Allocating variables for OzoneLunaMod 
    allocate(this%o3coefjmaxsha_patch(begp:endp))  ; this%o3coefjmaxsha_patch(:) = nan
    allocate(this%o3coefjmaxsun_patch(begp:endp))  ; this%o3coefjmaxsun_patch(:) = nan
    

  end subroutine InitAllocateLuna

 !-----------------------------------------------------------------------
  subroutine InitHistoryLuna(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize ozone history variables
    !
    ! !USES:
    use histFileMod  , only : hist_addfld1d
    !
    ! !ARGUMENTS:
    class(ozone_luna_type), intent(inout) :: this
    type(bounds_type), intent(in)    :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    
    character(len=*), parameter :: subname = 'InitHistoryLuna'
    !-----------------------------------------------------------------------
    
    begp = bounds%begp
    endp = bounds%endp

    call this%InitHistory(bounds)
    
    ! JmaxO3/Jmax0
    ! write description
    this%o3coefjmaxsun_patch(begp:endp) = spval
    call hist_addfld1d (fname='O3COEFJMAXSUN', units='1', &
         avgflag='A', long_name='LUNA ozone coefficient for Jmax for sunlit leaves', &
         ptr_patch=this%o3coefjmaxsun_patch)

    this%o3coefjmaxsha_patch(begp:endp) = spval
    call hist_addfld1d (fname='O3COEFJMAXSHA', units='1', &
         avgflag='A', long_name='LUNA ozone coefficient for Jmax for shaded leaves', &
         ptr_patch=this%o3coefjmaxsha_patch)

  end subroutine InitHistoryLuna
 
!-----------------------------------------------------------------------
  subroutine InitColdLuna(this, bounds)
    !
    ! !DESCRIPTION:
    ! Perform cold-start initialization for ozone
    !
    ! !ARGUMENTS:
    class(ozone_luna_type), intent(inout) :: this
    type(bounds_type), intent(in)    :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    
    character(len=*), parameter :: subname = 'InitColdLuna'
    !-----------------------------------------------------------------------

    begp = bounds%begp
    endp = bounds%endp

    call this%InitCold(bounds)
    
    this%o3coefjmaxsha_patch(begp:endp) = 0._r8
    this%o3coefjmaxsun_patch(begp:endp) = 0._r8
    
  end subroutine InitColdLuna

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
    class(ozone_luna_type) :: this
    type(bounds_type), intent(in)    :: bounds
    type(file_desc_t), intent(inout) :: ncid   ! netcdf id
    character(len=*) , intent(in)    :: flag   ! 'read', 'write' or 'define'
    !
    ! !LOCAL VARIABLES:
    logical :: readvar

    character(len=*), parameter :: subname = 'Restart'
    !-----------------------------------------------------------------------
    ! From ozone_type
    call restartvar(ncid=ncid, flag=flag, varname='o3_tlaiold', xtype=ncd_double, &
         dim1name='pft', &
         long_name='one-sided leaf area index, from previous timestep, for ozone calculations', units='', &
         readvar=readvar, interpinic_flag='interp', data=this%tlai_old_patch)
    
    call restartvar(ncid=ncid, flag=flag, varname='o3uptakesha', xtype=ncd_double, &
         dim1name='pft', &
         long_name='accumulated ozone uptake for shaded leaves', units='mmol m^-3', &
         readvar=readvar, interpinic_flag='interp', data=this%o3uptakesha_patch)

    call restartvar(ncid=ncid, flag=flag, varname='o3uptakesun', xtype=ncd_double, &
         dim1name='pft', &
         long_name='accumulated ozone uptake for sunlit leaves', units='mmol m^-3', &
         readvar=readvar, interpinic_flag='interp', data=this%o3uptakesun_patch)

    call restartvar(ncid=ncid, flag=flag, varname='o3coefjmaxsun', xtype=ncd_double, &
         dim1name='pft', &
         long_name='ozone coefficient for max. electron transport rate for sunlit leaves', units='unitless', &
         readvar=readvar, interpinic_flag='interp', data=this%o3coefjmaxsun_patch)

    call restartvar(ncid=ncid, flag=flag, varname='o3coefjmaxsha', xtype=ncd_double, &
         dim1name='pft', &
         long_name='ozone coefficient for max. electron transport rate for shaded leaves', units='unitless', &
         readvar=readvar, interpinic_flag='interp', data=this%o3coefjmaxsha_patch)

  end subroutine Restart

  ! ========================================================================
  ! Science routines
  ! ========================================================================

  !-----------------------------------------------------------------------
  subroutine CalcOzoneStress(this, bounds, num_exposedvegp, filter_exposedvegp, &
       forc_po3, forc_pbot, forc_th, rssun, rssha, rb, ram, tlai)
    !
    ! !DESCRIPTION:
    ! Calculate ozone uptake for each timestep.
    ! For consistency with the calls from CanopyFluxes, ozone coefficiens g/v are set to 1, respectively.
    !
    ! !USES:
    use PatchType            , only : patch
    !
    ! !ARGUMENTS:
    class(ozone_luna_type) , intent(inout) :: this
    type(bounds_type)      , intent(in)    :: bounds
    integer  , intent(in) :: num_exposedvegp           ! number of points in filter_exposedvegp
    integer  , intent(in) :: filter_exposedvegp(:)     ! patch filter for non-snow-covered veg
    real(r8) , intent(in) :: forc_po3( bounds%begc: )  ! ozone concentration in partial pressure (Pa)
    real(r8) , intent(in) :: forc_pbot( bounds%begc: ) ! atmospheric pressure (Pa)
    real(r8) , intent(in) :: forc_th( bounds%begc: )   ! atmospheric potential temperature (K)
    real(r8) , intent(in) :: rssun( bounds%begp: )     ! leaf stomatal resistance, sunlit leaves (s/m)
    real(r8) , intent(in) :: rssha( bounds%begp: )     ! leaf stomatal resistance, shaded leaves (s/m)
    real(r8) , intent(in) :: rb( bounds%begp: )        ! boundary layer resistance (s/m)
    real(r8) , intent(in) :: ram( bounds%begp: )       ! aerodynamical resistance (s/m)
    real(r8) , intent(in) :: tlai( bounds%begp: )      ! one-sided leaf area index, no burying by snow
    ! !LOCAL VARIABLES:
    integer  :: fp             ! filter index
    integer  :: p              ! patch index
    integer  :: c              ! column index

    real(r8) :: forc_ozone     ! ozone partial pressure (mol/mol)

    character(len=*), parameter :: subname = 'CalcOzoneStress'
    !-----------------------------------------------------------------------
    
    ! Enforce expected array sizes (mainly so that a debug-mode threaded test with
    ! ozone-off can pick up problems with the call to this routine)
    SHR_ASSERT_ALL_FL((ubound(forc_po3) == (/bounds%endc/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(forc_pbot) == (/bounds%endc/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(forc_th) == (/bounds%endc/)), sourcefile, __LINE__)
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

       ! Need to call to CalcOzoneUptake to update accumulated ozone over a days period.
       ! Convert from Pa to mol/mol
       forc_ozone = forc_po3(c) / forc_pbot(c)

       ! Ozone uptake for shaded leaves
       call this%CalcOzoneUptakeOnePoint( &
            forc_ozone=forc_ozone, forc_pbot=forc_pbot(c), forc_th=forc_th(c), &
            rs=rssha(p), rb=rb(p), ram=ram(p), &
            tlai=tlai(p), tlai_old=tlai_old(p), pft_type=patch%itype(p), &
            o3uptake=o3uptakesha(p))

       ! Ozone uptake for sunlit leaves
       call this%CalcOzoneUptakeOnePoint( &
            forc_ozone=forc_ozone, forc_pbot=forc_pbot(c), forc_th=forc_th(c), &
            rs=rssun(p), rb=rb(p), ram=ram(p), &
            tlai=tlai(p), tlai_old=tlai_old(p), pft_type=patch%itype(p), &
            o3uptake=o3uptakesun(p))

       tlai_old(p) = tlai(p)

    end do

    end associate

    ! Explicitly set outputs to 1. This isn't really needed, because they should still be
    ! at 1 from cold-start initialization, but do this for clarity here.

    this%o3coefvsha_patch(bounds%begp:bounds%endp) = 1._r8
    this%o3coefvsun_patch(bounds%begp:bounds%endp) = 1._r8
    this%o3coefgsha_patch(bounds%begp:bounds%endp) = 1._r8
    this%o3coefgsun_patch(bounds%begp:bounds%endp) = 1._r8

  end subroutine CalcOzoneStress

!-----------------------------------------------------------------------  
  subroutine Acc24_OzoneStressOnePoint_Luna(pft_type, o3uptake, o3coefjmax)

    ! !USES:
    use pftconMod            , only : pftcon
    !
    ! ! ARGUMENTS:
    integer  , intent(in)    :: pft_type    ! vegetation type, for indexing into pftvarcon arrays
    real(r8) , intent(in)    :: o3uptake    ! ozone entering the leaf
    real(r8) , intent(inout) :: o3coefjmax  ! ozone coefficient for max. electron transport rate
    !
    ! !LOCAL VARIABLES:
    real(r8) :: jmaxInt        ! intercept for max. electron transport rate
    real(r8) :: jmaxSlope      ! slope for max. electron transport rate

    if ( o3uptake == 0._r8 ) then
       ! No o3 damage if no o3 uptake
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


  end subroutine Acc24_OzoneStressOnePoint_Luna

!-----------------------------------------------------------------------  
  subroutine Acc24_OzoneStress_Luna(this, bounds, num_exposedvegp, filter_exposedvegp)
    !
    ! !DESCRIPTION:
    ! Calculate accumulated 24h ozone stress for luna.
    !
    ! !USES:
    use PatchType            , only : patch
    use pftconMod            , only : pftcon
    !
    ! !ARGUMENTS:
    class(ozone_luna_type) , intent(inout) :: this
    type(bounds_type)      , intent(in)    :: bounds
    integer  , intent(in)    :: num_exposedvegp        ! number of points in filter_exposedvegp
    integer  , intent(in)    :: filter_exposedvegp(:)  ! patch filter for non-snow-covered veg
    !
    ! !LOCAL VARIABLES:
    integer  :: fp             ! filter index
    integer  :: p              ! patch index
    integer  :: c              ! column index

    character(len=*), parameter :: subname = 'Acc24_OzoneStress_Luna'
    !-----------------------------------------------------------------------

    associate( &
         o3uptakesha => this%o3uptakesha_patch                , & ! Output: [real(r8) (:)] ozone dose
         o3uptakesun => this%o3uptakesun_patch                , & ! Output: [real(r8) (:)] ozone dose
         o3coefjmaxsha => this%o3coefjmaxsha_patch            , & ! Output: [real(r8) (:)] ozone coef jmax sha
         o3coefjmaxsun => this%o3coefjmaxsun_patch              & ! Output: [real(r8) (:)] ozone coef jmax sun
         )
      
      do fp = 1, num_exposedvegp
         p = filter_exposedvegp(fp)
         c = patch%column(p)

         ! Ozone damage for shaded leaves
         call this%Acc24_OzoneStressOnePoint_Luna( &
              pft_type=patch%itype(p), o3uptake=o3uptakesun(p), &
              o3coefjmax=o3coefjmaxsha(p))

         ! Ozone damage for sunlit leaves
         call this%Acc24_OzoneStressOnePoint_Luna( &
              pft_type=patch%itype(p), o3uptake=o3uptakesun(p), &
              o3coefjmax=o3coefjmaxsun(p))


      end do
    end associate

  end subroutine Acc24_OzoneStress_Luna

end module OzoneLunaMod
