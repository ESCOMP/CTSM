set (pfunit_sources
  test_ComputeSeedAmounts.pf)

create_pFUnit_test(CNVegComputeSeed test_CNVegComputeSeed_exe
  "${pfunit_sources}" "")

target_link_libraries(test_CNVegComputeSeed_exe clm csm_share esmf_wrf_timemgr)