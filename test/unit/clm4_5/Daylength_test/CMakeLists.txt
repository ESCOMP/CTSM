set(pf_sources test_daylength.pf)

set(all_sources "${share_sources}" "${clm_sources}")

set(sources_needed
  DaylengthMod.F90 clmtype.F90
  shr_kind_mod.F90 shr_infnan_mod.F90 shr_const_mod.F90
  decompMod_boundsTypeDecl.F90)

extract_sources("${sources_needed}" "${all_sources}" test_sources)

add_pFUnit_executable(test_daylength_exe
  "${pf_sources}"
  ${CMAKE_CURRENT_BINARY_DIR}
  "${test_sources}")

declare_generated_dependencies(test_daylength_exe "${share_genf90_sources}")

add_test(daylength test_daylength_exe)

define_pFUnit_failure(daylength)
