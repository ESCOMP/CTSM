module WaterBalanceType

  ! This type contains variables used in water balance checks and adjustments done based
  ! on changes in water states from one point to another. We'll need separate instances
  ! for each isotope/tracer if we want to do balance checks and related adjustments on
  ! the isotopes/tracers.
  type, public :: waterbalance_type
     real(r8), pointer :: h2osno_old_col         (:)   ! col snow mass for previous time step (kg/m2) (new)
     real(r8), pointer :: liq1_grc               (:)   ! grc initial gridcell total h2o liq content
     real(r8), pointer :: liq2_grc               (:)   ! grc post land cover change total liq content
     real(r8), pointer :: ice1_grc               (:)   ! grc initial gridcell total h2o ice content
     real(r8), pointer :: ice2_grc               (:)   ! grc post land cover change total ice content

     real(r8), pointer :: begwb_col              (:)   ! water mass begining of the time step
     real(r8), pointer :: endwb_col              (:)   ! water mass end of the time step
     real(r8), pointer :: errh2o_patch           (:)   ! water conservation error (mm H2O)
     real(r8), pointer :: errh2o_col             (:)   ! water conservation error (mm H2O)
     real(r8), pointer :: errh2osno_col          (:)   ! snow water conservation error(mm H2O)
   contains
     procedure :: Init
     procedure :: Restart
  end type waterbalance_type

contains

  ! Standard infrastructure routines here

end module WaterBalanceType
