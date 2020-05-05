set (pfunit_sources
  test_daylength.pf
  test_compute_max_daylength.pf)

create_pFUnit_test(Daylength test_Daylength_exe
  "${pfunit_sources}" "")

target_link_libraries(test_Daylength_exe clm csm_share esmf_wrf_timemgr)