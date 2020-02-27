module ctsm_LilacCouplingFieldIndices

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Defines all possible coupling field indices for coupling between atmosphere and land
  !
  ! !USES:
  use lilac_constants, only : field_index_unset

  implicit none
  private
  !
  ! !PUBLIC DATA:

  ! ------------------------------------------------------------------------
  ! These are the fields that can be passed from atm -> lnd. The host atmosphere model
  ! will refer to these indices when setting fields.
  ! ------------------------------------------------------------------------

  integer, public :: lilac_a2l_Sa_landfrac = field_index_unset
  integer, public :: lilac_a2l_Sa_z = field_index_unset
  integer, public :: lilac_a2l_Sa_topo = field_index_unset
  integer, public :: lilac_a2l_Sa_u = field_index_unset
  integer, public :: lilac_a2l_Sa_v = field_index_unset
  integer, public :: lilac_a2l_Sa_ptem = field_index_unset
  integer, public :: lilac_a2l_Sa_pbot = field_index_unset
  integer, public :: lilac_a2l_Sa_tbot = field_index_unset
  integer, public :: lilac_a2l_Sa_shum = field_index_unset
  integer, public :: lilac_a2l_Faxa_lwdn = field_index_unset
  integer, public :: lilac_a2l_Faxa_rainc = field_index_unset
  integer, public :: lilac_a2l_Faxa_rainl = field_index_unset
  integer, public :: lilac_a2l_Faxa_snowc = field_index_unset
  integer, public :: lilac_a2l_Faxa_snowl = field_index_unset
  integer, public :: lilac_a2l_Faxa_swndr = field_index_unset
  integer, public :: lilac_a2l_Faxa_swvdr = field_index_unset
  integer, public :: lilac_a2l_Faxa_swndf = field_index_unset
  integer, public :: lilac_a2l_Faxa_swvdf = field_index_unset

  integer, public :: lilac_a2l_Faxa_bcphidry = field_index_unset
  integer, public :: lilac_a2l_Faxa_bcphodry = field_index_unset
  integer, public :: lilac_a2l_Faxa_bcphiwet = field_index_unset
  integer, public :: lilac_a2l_Faxa_ocphidry = field_index_unset
  integer, public :: lilac_a2l_Faxa_ocphodry = field_index_unset
  integer, public :: lilac_a2l_Faxa_ocphiwet = field_index_unset
  integer, public :: lilac_a2l_Faxa_dstwet1 = field_index_unset
  integer, public :: lilac_a2l_Faxa_dstdry1 = field_index_unset
  integer, public :: lilac_a2l_Faxa_dstwet2 = field_index_unset
  integer, public :: lilac_a2l_Faxa_dstdry2 = field_index_unset
  integer, public :: lilac_a2l_Faxa_dstwet3 = field_index_unset
  integer, public :: lilac_a2l_Faxa_dstdry3 = field_index_unset
  integer, public :: lilac_a2l_Faxa_dstwet4 = field_index_unset
  integer, public :: lilac_a2l_Faxa_dstdry4 = field_index_unset

  ! ------------------------------------------------------------------------
  ! These are the fields that can be passed from lnd -> atm. The host atmosphere model
  ! will refer to these indices when retrieving fields.
  ! ------------------------------------------------------------------------

  integer, public :: lilac_l2a_Sl_t = field_index_unset
  integer, public :: lilac_l2a_Sl_tref = field_index_unset
  integer, public :: lilac_l2a_Sl_qref = field_index_unset
  integer, public :: lilac_l2a_Sl_avsdr = field_index_unset
  integer, public :: lilac_l2a_Sl_anidr = field_index_unset
  integer, public :: lilac_l2a_Sl_avsdf = field_index_unset
  integer, public :: lilac_l2a_Sl_anidf = field_index_unset
  integer, public :: lilac_l2a_Sl_snowh = field_index_unset
  integer, public :: lilac_l2a_Sl_u10 = field_index_unset
  integer, public :: lilac_l2a_Sl_fv = field_index_unset
  integer, public :: lilac_l2a_Sl_ram1 = field_index_unset
  integer, public :: lilac_l2a_Sl_z0m = field_index_unset
  integer, public :: lilac_l2a_Fall_taux = field_index_unset
  integer, public :: lilac_l2a_Fall_tauy = field_index_unset
  integer, public :: lilac_l2a_Fall_lat = field_index_unset
  integer, public :: lilac_l2a_Fall_sen = field_index_unset
  integer, public :: lilac_l2a_Fall_lwup = field_index_unset
  integer, public :: lilac_l2a_Fall_evap = field_index_unset
  integer, public :: lilac_l2a_Fall_swnet = field_index_unset
  integer, public :: lilac_l2a_Fall_flxdst1 = field_index_unset
  integer, public :: lilac_l2a_Fall_flxdst2 = field_index_unset
  integer, public :: lilac_l2a_Fall_flxdst3 = field_index_unset
  integer, public :: lilac_l2a_Fall_flxdst4 = field_index_unset

end module ctsm_LilacCouplingFieldIndices
