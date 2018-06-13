module WaterStateType

  ! This type contains water state variables for which we need a separate instance for
  ! each isotope or water tracer.
  type, public :: waterstate_type
     real(r8), pointer :: h2osno_col             (:)   ! col snow water (mm H2O)
     real(r8), pointer :: h2osoi_liq_col         (:,:) ! col liquid water (kg/m2) (new) (-nlevsno+1:nlevgrnd)    
     real(r8), pointer :: h2osoi_ice_col         (:,:) ! col ice lens (kg/m2) (new) (-nlevsno+1:nlevgrnd)    
     real(r8), pointer :: h2ocan_patch           (:)   ! patch canopy water (mm H2O)  NOTE: Will eventually be moved to diagnostic (see https://github.com/ESCOMP/ctsm/issues/199)
     real(r8), pointer :: h2osfc_col             (:)   ! col surface water (mm H2O)
     real(r8), pointer :: snocan_patch           (:)   ! patch canopy snow water (mm H2O)
     real(r8), pointer :: liqcan_patch           (:)   ! patch canopy liquid water (mm H2O)
   contains
     procedure :: Init
     procedure :: Restart
  end type waterstate_type

contains

  ! Standard infrastructure routines here

end module WaterStateType
