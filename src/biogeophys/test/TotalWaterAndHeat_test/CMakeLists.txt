set (pfunit_sources
  test_total_water_and_heat.pf)

add_pfunit_ctest(total_water_and_heat
  TEST_SOURCES "${pfunit_sources}"
  LINK_LIBRARIES clm csm_share esmf
  EXTRA_FINALIZE unittest_finalize_esmf
  EXTRA_USE unittestInitializeAndFinalize)
