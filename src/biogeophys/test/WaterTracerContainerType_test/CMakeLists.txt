set (pfunit_sources
  test_water_tracer_container.pf)

# extra sources used for this test, which are not .pf files
# (currently none)
set (extra_sources
  )

create_pFUnit_test(water_tracer_container test_water_tracer_container_exe
  "${pfunit_sources}" "${extra_sources}")

target_link_libraries(test_water_tracer_container_exe clm csm_share esmf_wrf_timemgr)