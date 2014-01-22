set(all_sources "${share_sources}" "${clm_sources}")

set(sources_needed
  DaylengthMod.F90 clmtype.F90
  shr_kind_mod.F90 shr_infnan_mod.F90 shr_const_mod.F90
  decompMod_boundsTypeDecl.F90)

extract_sources("${sources_needed}" "${all_sources}" test_sources)

create_pFUnit_test(daylength test_daylength_exe
  "test_daylength.pf" "${test_sources}")
