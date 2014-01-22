set(all_sources "${share_sources}" "${clm_sources}")

set(sources_needed
  dynTimeInfoMod.F90 clm_varctl.F90 
  shr_assert_mod.F90 shr_kind_mod.F90 shr_infnan_mod.F90 shr_log_mod.F90 shr_sys_mod_mock.F90)

extract_sources("${sources_needed}" "${all_sources}" test_sources)

create_pFUnit_test(dynTimeInfo test_dynTimeInfo_exe
  "test_dynTimeInfo.pf" "${test_sources}")

