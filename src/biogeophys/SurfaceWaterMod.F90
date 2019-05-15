module SurfaceWaterMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Routines for handling surface water (h2osfc) and related terms
  !
  ! !USES:
#include "shr_assert.h"
  use shr_kind_mod            , only : r8 => shr_kind_r8
  use shr_log_mod             , only : errMsg => shr_log_errMsg
  use shr_const_mod           , only : shr_const_pi
  use shr_spfn_mod            , only : erf => shr_spfn_erf
  use landunit_varcon         , only : istsoil, istcrop
  use decompMod               , only : bounds_type
  use ColumnType              , only : col
  use LandunitType            , only : lun
  use WaterStateBulkType      , only : waterstatebulk_type
  use WaterDiagnosticBulkType , only : waterdiagnosticbulk_type

  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: FracH2oSfc           ! Determine fraction of land surfaces which are submerged  

contains

  !-----------------------------------------------------------------------
  subroutine FracH2oSfc(bounds, num_h2osfc, filter_h2osfc, &
       waterstatebulk_inst, waterdiagnosticbulk_inst, no_update)
    !
    ! !DESCRIPTION:
    ! Determine fraction of land surfaces which are submerged  
    ! based on surface microtopography and surface water storage.
    !
    ! !ARGUMENTS:
    type(bounds_type)     , intent(in)           :: bounds           
    integer               , intent(in)           :: num_h2osfc       ! number of column points in column filter
    integer               , intent(in)           :: filter_h2osfc(:) ! column filter 
    type(waterstatebulk_type) , intent(inout)        :: waterstatebulk_inst
    type(waterdiagnosticbulk_type) , intent(inout)        :: waterdiagnosticbulk_inst
    integer               , intent(in), optional :: no_update        ! flag to make calculation w/o updating variables
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

    do f = 1, num_h2osfc
       c = filter_h2osfc(f)
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


          if (.not. present(no_update)) then

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

          endif ! end of no_update construct

       else !if landunit not istsoil/istcrop, set frac_h2osfc to zero

          frac_h2osfc(c) = 0._r8

       endif

    end do

    end associate

  end subroutine FracH2oSfc

end module SurfaceWaterMod
