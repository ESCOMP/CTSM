set (pfunit_sources
  test_FATESFireFactory.pf
)

add_pfunit_ctest(FATESFireFActory
  TEST_SOURCES "${pfunit_sources}"
  LINK_LIBRARIES clm csm_share esmf)
