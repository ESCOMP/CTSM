module WaterDiagnosticType

  ! This type contains water diagnostic variables for which we need a separate instance
  ! for each isotope or water tracer.
  type, public :: waterdiag_type
     real(r8), pointer :: snowice_col            (:)   ! col average snow ice lens
     real(r8), pointer :: snowliq_col            (:)   ! col average snow liquid water
     real(r8), pointer :: h2osoi_liqice_10cm_col (:)   ! col liquid water + ice lens in top 10cm of soil (kg/m2)
     real(r8), pointer :: h2osoi_vol_col         (:,:) ! col volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]  (nlevgrnd)

     real(r8), pointer :: q_ref2m_patch          (:)   ! patch 2 m height surface specific humidity (kg/kg)
     real(r8), pointer :: qg_col                 (:)   ! col ground specific humidity [kg/kg]
     real(r8), pointer :: tws_grc                (:)   ! grc total water storage (mm H2O)  NOTE: Similar to wt_col in the old code

     ! NOTE(wjs, 2018-05-24) I'm not sure whether we need a separate instance for each
     ! tracer for these, but they look similar to qg, which had a separate water isotope
     ! instance in the old water isotope branch.
     real(r8), pointer :: qg_snow_col            (:)   ! col ground specific humidity [kg/kg]
     real(r8), pointer :: qg_soil_col            (:)   ! col ground specific humidity [kg/kg]
     real(r8), pointer :: qg_h2osfc_col          (:)   ! col ground specific humidity [kg/kg]
     real(r8), pointer :: qaf_lun                (:)   ! lun urban canopy air specific humidity (kg/kg)
   contains
     procedure :: Init
     procedure :: Restart
  end type waterdiag_type

contains

  ! Standard infrastructure routines here

end module WaterDiagnosticType
