set (pfunit_sources
  test_SpeciesNonIsotope.pf
  test_SpeciesIsotope.pf)

create_pFUnit_test(Species test_Species_exe
  "${pfunit_sources}" "")

target_link_libraries(test_Species_exe clm csm_share)