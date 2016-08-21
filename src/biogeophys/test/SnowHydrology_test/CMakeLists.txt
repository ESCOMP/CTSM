set (pfunit_sources
  test_SnowHydrology_initSnowLayers.pf
  test_SnowHydrology_newSnowBulkDensity.pf)

create_pFUnit_test(SnowHydrology test_SnowHydrology_exe
  "${pfunit_sources}" "")

target_link_libraries(test_SnowHydrology_exe clm csm_share esmf_wrf_timemgr)
