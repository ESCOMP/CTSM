set (pfunit_sources
  test_calc_tracer_from_bulk_fixed_ratio.pf
  test_calc_tracer_from_bulk.pf
  test_compare_bulk_to_tracer.pf)

add_pfunit_ctest(water_tracer_utils
  TEST_SOURCES "${pfunit_sources}"
  LINK_LIBRARIES clm csm_share esmf
  EXTRA_FINALIZE unittest_finalize_esmf
  EXTRA_USE unittestInitializeAndFinalize)
