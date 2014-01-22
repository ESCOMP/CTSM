set(all_sources "${share_sources}" "${clm_sources}")

set (pfunit_sources
  test_dynVarTimeUninterp.pf
  test_dynVarTimeInterp.pf)

set (sources_needed
  dynVarMod.F90 dynVarTimeInterpMod.F90 dynVarTimeUninterpMod.F90
  dynTimeInfoMod.F90 dynFileMod_mock.F90
  surfrdUtilsMod.F90 clm_varctl.F90 clm_time_manager_mock.F90 ncdio_pio_mock.F90 ncdio_var.F90 spmdMod_mock.F90
  shr_assert_mod.F90 shr_kind_mod.F90 shr_infnan_mod.F90 shr_log_mod.F90 shr_sys_mod_mock.F90)

extract_sources("${sources_needed}" "${all_sources}" test_sources)

# extra sources used for this test, which are not .pf files, but also aren't obtained from the all_sources list
set (extra_sources
  test_dynVarShared.F90)

create_pfUnit_test(dynVar test_dynVar_exe
  "${pfunit_sources}" "${test_sources};${extra_sources}")

# This test tests some modules generated with genf90, so we used declare_generated_dependencies to specify that
declare_generated_dependencies(test_dynVar_exe "${clm_genf90_sources}")
