set(pfunit_sources
  test_accumul.pf)

create_pFUnit_test(accumul test_accumul_exe
  "${pfunit_sources}" "")

target_link_libraries(test_accumul_exe clm csm_share esmf_wrf_timemgr)