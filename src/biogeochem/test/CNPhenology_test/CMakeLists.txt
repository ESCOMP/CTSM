set (pfunit_sources
  test_CNPhenology.pf)

create_pFUnit_test(CNPhenology test_CNPhenology_exe
  "${pfunit_sources}" "")

target_link_libraries(test_CNPhenology_exe clm csm_share esmf_wrf_timemgr)
