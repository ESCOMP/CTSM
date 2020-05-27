set(pfunit_sources
  test_dyn_cons_biogeophys.pf)

create_pFUnit_test(dynConsBiogeophys test_dynConsBiogeophys_exe
  "${pfunit_sources}" "")

target_link_libraries(test_dynConsBiogeophys_exe clm csm_share esmf_wrf_timemgr)

