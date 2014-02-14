set(all_sources "${share_sources}" "${clm_sources}")

set(pfunit_sources
  test_update_landunit_weights_one_gcell.pf
  test_update_landunit_weights.pf)

set(sources_needed
  dynLandunitAreaMod.F90 clmtype.F90 clm_varctl.F90 clm_varcon.F90 clm_varpar.F90
  decompMod_boundsTypeDecl.F90
  shr_assert_mod.F90 shr_kind_mod.F90 shr_const_mod.F90 shr_infnan_mod.F90 shr_log_mod.F90
  shr_sys_mod_mock.F90)

extract_sources("${sources_needed}" "${all_sources}" test_sources)

create_pFUnit_test(dynLandunitArea test_dynLandunitArea_exe
  "${pfunit_sources}" "${test_sources}")

