set (pfunit_sources
  test_CNFireFactory.pf
)

add_pfunit_ctest(CNFireFActory
  TEST_SOURCES "${pfunit_sources}"
  LINK_LIBRARIES clm csm_share esmf)
