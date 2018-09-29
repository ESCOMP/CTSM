set (pfunit_sources
  test_calc_tracer_from_bulk_fixed_ratio.pf
  test_calc_tracer_from_bulk.pf
  test_compare_bulk_to_tracer.pf)

# extra sources used for this test, which are not .pf files
# (currently none)
set (extra_sources
  )

create_pFUnit_test(water_tracer_utils test_water_tracer_utils_exe
  "${pfunit_sources}" "${extra_sources}")

target_link_libraries(test_water_tracer_utils_exe clm csm_share esmf_wrf_timemgr)