module SurfaceWaterMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Routines for handling surface water (h2osfc) and related terms
  !
  ! !USES:
#include "shr_assert.h"
  use shr_kind_mod                , only : r8 => shr_kind_r8
  use shr_log_mod                 , only : errMsg => shr_log_errMsg
  use shr_const_mod               , only : shr_const_pi
  use shr_spfn_mod                , only : erf => shr_spfn_erf
  use clm_varcon                  , only : denh2o, denice, roverg, wimp, tfrz, pc, mu, rpi
  use clm_time_manager            , only : get_step_size
  use column_varcon               , only : icol_roof, icol_road_imperv, icol_sunwall, icol_shadewall, icol_road_perv
  use landunit_varcon             , only : istsoil, istcrop
  use decompMod                   , only : bounds_type
  use ColumnType                  , only : col
  use LandunitType                , only : lun
  use NumericsMod                 , only : truncate_small_values
  use InfiltrationExcessRunoffMod , only : infiltration_excess_runoff_type
  use EnergyFluxType              , only : energyflux_type
  use SoilHydrologyType           , only : soilhydrology_type
  use WaterFluxBulkType           , only : waterfluxbulk_type
  use WaterStateBulkType          , only : waterstatebulk_type
  use WaterDiagnosticBulkType     , only : waterdiagnosticbulk_type

  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: FracH2oSfc           ! Determine fraction of land surfaces which are submerged  
  public :: UpdateH2osfc         ! Calculate fluxes out of h2osfc and update the h2osfc state

  ! !PRIVATE MEMBER FUNCTIONS:
  private :: QflxH2osfcSurf      ! Compute qflx_h2osfc_surf
  private :: QflxH2osfcDrain     ! Compute qflx_h2osfc_drain

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

contains

  !-----------------------------------------------------------------------
  subroutine FracH2oSfc(bounds, num_nolakec, filter_nolakec, &
       waterstatebulk_inst, waterdiagnosticbulk_inst)
    !
    ! !DESCRIPTION:
    ! Determine fraction of land surfaces which are submerged  
    ! based on surface microtopography and surface water storage.
    !
    ! !ARGUMENTS:
    type(bounds_type)     , intent(in)           :: bounds           
    integer               , intent(in)           :: num_nolakec       ! number of points in nolakec filter
    integer               , intent(in)           :: filter_nolakec(:) ! column filter for non-lake points 
    type(waterstatebulk_type) , intent(inout)        :: waterstatebulk_inst
    type(waterdiagnosticbulk_type) , intent(inout)        :: waterdiagnosticbulk_inst
    !
    ! !LOCAL VARIABLES:
    integer :: c,f,l          ! indices
    real(r8):: d,fd,dfdd      ! temporary variable for frac_h2o iteration
    real(r8):: sigma          ! microtopography pdf sigma in mm
    real(r8):: min_h2osfc
    !-----------------------------------------------------------------------

    associate(                                              & 
         micro_sigma  => col%micro_sigma                  , & ! Input:  [real(r8) (:)   ] microtopography pdf sigma (m)                     
         
         h2osno       => waterstatebulk_inst%h2osno_col       , & ! Input:  [real(r8) (:)   ] snow water (mm H2O)                               
         
         h2osoi_liq   => waterstatebulk_inst%h2osoi_liq_col   , & ! Output: [real(r8) (:,:) ] liquid water (col,lyr) [kg/m2]                  
         h2osfc       => waterstatebulk_inst%h2osfc_col       , & ! Output: [real(r8) (:)   ] surface water (mm)                                
         frac_sno     => waterdiagnosticbulk_inst%frac_sno_col     , & ! Output: [real(r8) (:)   ] fraction of ground covered by snow (0 to 1)       
         frac_sno_eff => waterdiagnosticbulk_inst%frac_sno_eff_col , & ! Output: [real(r8) (:)   ] eff. fraction of ground covered by snow (0 to 1)  
         frac_h2osfc  => waterdiagnosticbulk_inst%frac_h2osfc_col,   & ! Output: [real(r8) (:)   ] col fractional area with surface water greater than zero 
         frac_h2osfc_nosnow  => waterdiagnosticbulk_inst%frac_h2osfc_nosnow_col    & ! Output: [real(r8) (:)   ] col fractional area with surface water greater than zero (if no snow present)
         )

    ! arbitrary lower limit on h2osfc for safer numerics...
    min_h2osfc=1.e-8_r8

    do f = 1, num_nolakec
       c = filter_nolakec(f)
       l = col%landunit(c)

       ! h2osfc only calculated for soil vegetated land units
       if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then

          !  Use newton-raphson method to iteratively determine frac_h2osfc
          !  based on amount of surface water storage (h2osfc) and 
          !  microtopography variability (micro_sigma)

          if (h2osfc(c) > min_h2osfc) then
             ! a cutoff is needed for numerical reasons...(nonconvergence after 5 iterations)
             d=0.0

             sigma=1.0e3 * micro_sigma(c) ! convert to mm
             do l=1,10
                fd = 0.5*d*(1.0_r8+erf(d/(sigma*sqrt(2.0)))) &
                     +sigma/sqrt(2.0*shr_const_pi)*exp(-d**2/(2.0*sigma**2)) &
                     -h2osfc(c)
                dfdd = 0.5*(1.0_r8+erf(d/(sigma*sqrt(2.0))))

                d = d - fd/dfdd
             enddo
             !--  update the submerged areal fraction using the new d value
             frac_h2osfc(c) = 0.5*(1.0_r8+erf(d/(sigma*sqrt(2.0))))

          else
             frac_h2osfc(c) = 0._r8
             h2osoi_liq(c,1) = h2osoi_liq(c,1) + h2osfc(c)
             h2osfc(c)=0._r8
          endif

          frac_h2osfc_nosnow(c) = frac_h2osfc(c)

          ! adjust fh2o, fsno when sum is greater than zero
          if (frac_sno(c) > (1._r8 - frac_h2osfc(c)) .and. h2osno(c) > 0) then

             if (frac_h2osfc(c) > 0.01_r8) then             
                frac_h2osfc(c) = max(1.0_r8 - frac_sno(c),0.01_r8)
                frac_sno(c) = 1.0_r8 - frac_h2osfc(c)
             else
                frac_sno(c) = 1.0_r8 - frac_h2osfc(c)
             endif
             frac_sno_eff(c)=frac_sno(c)

          endif

       else !if landunit not istsoil/istcrop, set frac_h2osfc to zero

          frac_h2osfc(c) = 0._r8

       endif

    end do

    end associate

  end subroutine FracH2oSfc

  !-----------------------------------------------------------------------
  subroutine UpdateH2osfc(bounds, num_hydrologyc, filter_hydrologyc, &
       infiltration_excess_runoff_inst, &
       energyflux_inst, soilhydrology_inst, &
       waterfluxbulk_inst, waterstatebulk_inst, waterdiagnosticbulk_inst)
    !
    ! !DESCRIPTION:
    ! Calculate fluxes out of h2osfc and update the h2osfc state
    !
    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)    :: bounds               
    integer                  , intent(in)    :: num_hydrologyc       ! number of column soil points in column filter
    integer                  , intent(in)    :: filter_hydrologyc(:) ! column filter for soil points
    type(infiltration_excess_runoff_type), intent(in) :: infiltration_excess_runoff_inst
    type(energyflux_type)    , intent(in)    :: energyflux_inst
    type(soilhydrology_type) , intent(in)    :: soilhydrology_inst
    type(waterfluxbulk_type)     , intent(inout) :: waterfluxbulk_inst
    type(waterstatebulk_type)    , intent(inout) :: waterstatebulk_inst
    type(waterdiagnosticbulk_type)    , intent(in) :: waterdiagnosticbulk_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: c,l,fc                                     ! indices
    real(r8) :: dtime                                      ! land model time step (sec)
    real(r8) :: h2osfc_partial(bounds%begc:bounds%endc)    ! partially-updated h2osfc
    !-----------------------------------------------------------------------

    associate(                                                        & 
         qinmax           =>    infiltration_excess_runoff_inst%qinmax_col , & ! Input:  [real(r8) (:)] maximum infiltration rate (mm H2O /s)
         
         frac_h2osfc      =>    waterdiagnosticbulk_inst%frac_h2osfc_col     , & ! Input:  [real(r8) (:)   ]  fraction of ground covered by surface water (0 to 1)
         frac_h2osfc_nosnow  => waterdiagnosticbulk_inst%frac_h2osfc_nosnow_col,    & ! Input: [real(r8) (:)   ] col fractional area with surface water greater than zero (if no snow present)
         h2osfc           =>    waterstatebulk_inst%h2osfc_col          , & ! Output: [real(r8) (:)   ]  surface water (mm)                                
         
         qflx_in_h2osfc   => waterfluxbulk_inst%qflx_in_h2osfc_col      , & ! Input:  [real(r8) (:)   ] total surface input to h2osfc
         qflx_h2osfc_surf =>    waterfluxbulk_inst%qflx_h2osfc_surf_col , & ! Output: [real(r8) (:)   ]  surface water runoff (mm H2O /s)                       
         qflx_h2osfc_drain => waterfluxbulk_inst%qflx_h2osfc_drain_col  , & ! Output: [real(r8) (:)   ]  bottom drainage from h2osfc (mm H2O /s)
         
         h2osfc_thresh    =>    soilhydrology_inst%h2osfc_thresh_col, & ! Input:  [real(r8) (:)   ]  level at which h2osfc "percolates"                
         h2osfcflag       =>    soilhydrology_inst%h2osfcflag         & ! Input:  integer
         )

    dtime = get_step_size()

    call QflxH2osfcSurf(bounds, num_hydrologyc, filter_hydrologyc, &
         h2osfcflag = h2osfcflag, &
         h2osfc = h2osfc(bounds%begc:bounds%endc), &
         h2osfc_thresh = h2osfc_thresh(bounds%begc:bounds%endc), &
         frac_h2osfc_nosnow = frac_h2osfc_nosnow(bounds%begc:bounds%endc), &
         topo_slope = col%topo_slope(bounds%begc:bounds%endc), &
         qflx_h2osfc_surf = qflx_h2osfc_surf(bounds%begc:bounds%endc))

    ! Update h2osfc prior to calculating bottom drainage from h2osfc.
    !
    ! This could be removed if we wanted to do a straight forward Euler, and/or set
    ! things up for a more sophisticated solution method. In the latter, rather than
    ! having h2osfc_partial, we'd have some more sophisticated state estimate here
    do fc = 1, num_hydrologyc
       c = filter_hydrologyc(fc)
       h2osfc_partial(c) = h2osfc(c) + (qflx_in_h2osfc(c) - qflx_h2osfc_surf(c)) * dtime
    end do

    call truncate_small_values(num_f = num_hydrologyc, filter_f = filter_hydrologyc, &
         lb = bounds%begc, ub = bounds%endc, &
         data_baseline = h2osfc(bounds%begc:bounds%endc), &
         data = h2osfc_partial(bounds%begc:bounds%endc))

    call QflxH2osfcDrain(bounds, num_hydrologyc, filter_hydrologyc, &
         h2osfcflag = h2osfcflag, &
         h2osfc = h2osfc_partial(bounds%begc:bounds%endc), &
         frac_h2osfc = frac_h2osfc(bounds%begc:bounds%endc), &
         qinmax = qinmax(bounds%begc:bounds%endc), &
         qflx_h2osfc_drain = qflx_h2osfc_drain(bounds%begc:bounds%endc))

    ! Update h2osfc based on fluxes
    do fc = 1, num_hydrologyc
       c = filter_hydrologyc(fc)
       h2osfc(c) = h2osfc_partial(c) - qflx_h2osfc_drain(c) * dtime
    end do

    ! Due to rounding errors, fluxes that should have brought h2osfc to exactly 0 may
    ! have instead left it slightly less than or slightly greater than 0. Correct for
    ! that here.
    call truncate_small_values(num_f = num_hydrologyc, filter_f = filter_hydrologyc, &
         lb = bounds%begc, ub = bounds%endc, &
         data_baseline = h2osfc_partial(bounds%begc:bounds%endc), &
         data = h2osfc(bounds%begc:bounds%endc))

    end associate

  end subroutine UpdateH2osfc

  !-----------------------------------------------------------------------
  subroutine QflxH2osfcSurf(bounds, num_hydrologyc, filter_hydrologyc, &
       h2osfcflag, h2osfc, h2osfc_thresh, frac_h2osfc_nosnow, topo_slope, &
       qflx_h2osfc_surf)
    !
    ! !DESCRIPTION:
    ! Compute qflx_h2osfc_surf
    !
    ! !ARGUMENTS:
    type(bounds_type) , intent(in)    :: bounds
    integer           , intent(in)    :: num_hydrologyc                     ! number of column soil points in column filter
    integer           , intent(in)    :: filter_hydrologyc(:)               ! column filter for soil points
    integer           , intent(in)    :: h2osfcflag                         ! true => surface water is active
    real(r8)          , intent(in)    :: h2osfc( bounds%begc: )             ! surface water (mm)
    real(r8)          , intent(in)    :: h2osfc_thresh( bounds%begc: )      ! level at which h2osfc "percolates"
    real(r8)          , intent(in)    :: frac_h2osfc_nosnow( bounds%begc: ) ! fractional area with surface water greater than zero (if no snow present)
    real(r8)          , intent(in)    :: topo_slope( bounds%begc: )         ! topographic slope
    real(r8)          , intent(inout) :: qflx_h2osfc_surf( bounds%begc: )   ! surface water runoff (mm H2O /s)
    !
    ! !LOCAL VARIABLES:
    integer  :: fc, c
    real(r8) :: dtime         ! land model time step (sec)
    real(r8) :: frac_infclust ! fraction of submerged area that is connected
    real(r8) :: k_wet         ! linear reservoir coefficient for h2osfc

    character(len=*), parameter :: subname = 'QflxH2osfcSurf'
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(h2osfc) == (/bounds%endc/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(h2osfc_thresh) == (/bounds%endc/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(frac_h2osfc_nosnow) == (/bounds%endc/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(topo_slope) == (/bounds%endc/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(qflx_h2osfc_surf) == (/bounds%endc/)), errMsg(sourcefile, __LINE__))

    dtime = get_step_size()

    do fc = 1, num_hydrologyc
       c = filter_hydrologyc(fc)

       if (h2osfcflag==1) then
          if (frac_h2osfc_nosnow(c) <= pc) then
             frac_infclust=0.0_r8
          else
             frac_infclust=(frac_h2osfc_nosnow(c)-pc)**mu
          endif
       endif

       ! limit runoff to value of storage above S(pc)
       if(h2osfc(c) > h2osfc_thresh(c) .and. h2osfcflag/=0) then
          ! spatially variable k_wet
          k_wet=1.0e-4_r8 * sin((rpi/180.) * topo_slope(c))
          qflx_h2osfc_surf(c) = k_wet * frac_infclust * (h2osfc(c) - h2osfc_thresh(c))

          qflx_h2osfc_surf(c)=min(qflx_h2osfc_surf(c),(h2osfc(c) - h2osfc_thresh(c))/dtime)
       else
          qflx_h2osfc_surf(c)= 0._r8
       endif

       ! cutoff lower limit
       if ( qflx_h2osfc_surf(c) < 1.0e-8) then
          qflx_h2osfc_surf(c) = 0._r8
       end if

    end do

  end subroutine QflxH2osfcSurf

  !-----------------------------------------------------------------------
  subroutine QflxH2osfcDrain(bounds, num_hydrologyc, filter_hydrologyc, &
       h2osfcflag, h2osfc, frac_h2osfc, qinmax, &
       qflx_h2osfc_drain)
    !
    ! !DESCRIPTION:
    ! Compute qflx_h2osfc_drain
    !
    ! Note that, if h2osfc is negative, then qflx_h2osfc_drain will be negative - acting
    ! to exactly restore h2osfc to 0.
    !
    ! !ARGUMENTS:
    type(bounds_type) , intent(in)    :: bounds
    integer           , intent(in)    :: num_hydrologyc                     ! number of column soil points in column filter
    integer           , intent(in)    :: filter_hydrologyc(:)               ! column filter for soil points
    integer           , intent(in)    :: h2osfcflag                         ! true => surface water is active
    real(r8)          , intent(in)    :: h2osfc( bounds%begc: )             ! surface water (mm)
    real(r8)          , intent(in)    :: frac_h2osfc( bounds%begc: )        ! fraction of ground covered by surface water (0 to 1)
    real(r8)          , intent(in)    :: qinmax( bounds%begc: )             ! maximum infiltration rate (mm H2O /s)
    real(r8)          , intent(inout) :: qflx_h2osfc_drain( bounds%begc: )  ! bottom drainage from h2osfc (mm H2O /s)
    !
    ! !LOCAL VARIABLES:
    integer :: fc, c
    real(r8) :: dtime         ! land model time step (sec)

    character(len=*), parameter :: subname = 'QflxH2osfcDrain'
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(h2osfc) == (/bounds%endc/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(frac_h2osfc) == (/bounds%endc/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(qinmax) == (/bounds%endc/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(qflx_h2osfc_drain) == (/bounds%endc/)), errMsg(sourcefile, __LINE__))

    dtime = get_step_size()

    do fc = 1, num_hydrologyc
       c = filter_hydrologyc(fc)

       if (h2osfc(c) < 0.0) then
          qflx_h2osfc_drain(c) = h2osfc(c)/dtime
       else
          qflx_h2osfc_drain(c)=min(frac_h2osfc(c)*qinmax(c),h2osfc(c)/dtime)
          if(h2osfcflag==0) then
             ! ensure no h2osfc
             qflx_h2osfc_drain(c)= max(0._r8,h2osfc(c)/dtime)
          end if
       end if
    end do

  end subroutine QflxH2osfcDrain

end module SurfaceWaterMod
