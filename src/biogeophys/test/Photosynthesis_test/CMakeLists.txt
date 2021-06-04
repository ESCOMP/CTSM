set (pfunit_sources
   test_Photosynthesis.pf)

create_pFUnit_test(Photosynthesis test_Photosynthesis_exe
  "${pfunit_sources}" "")

target_link_libraries(test_Photosynthesis_exe clm csm_share esmf_wrf_timemgr)
