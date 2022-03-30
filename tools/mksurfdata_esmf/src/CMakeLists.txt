cmake_minimum_required(VERSION 3.10)
project(mksurfdata Fortran)

set (ESMFMKFILE /glade/work/oehmke/ESMF/scalable_mesh_from_file/lib/libg/Linux.intel.64.mpt.default/esmf.mk)

file(STRINGS ${ESMFMKFILE} esmf_mk_text)
foreach(line ${esmf_mk_text})
  string(REGEX REPLACE "^[ ]+" "" line ${line}) # strip leading spaces
  if (line MATCHES "^ESMF_*")                   # process only line starting with ESMF_
    string(REGEX MATCH "^ESMF_[^=]+" esmf_name ${line})
    string(REPLACE "${esmf_name}=" "" emsf_value ${line})
    set(${esmf_name} "${emsf_value}")
  endif()
endforeach()
string(REPLACE "-I" "" ESMF_F90COMPILEPATHS ${ESMF_F90COMPILEPATHS})
string(REPLACE " " ";" ESMF_F90COMPILEPATHS ${ESMF_F90COMPILEPATHS})
message("ESMF_F90COMPILEPATHS:   ${ESMF_F90COMPILEPATHS}")

set(SRCFILES mkagfirepkmonthMod.F90
             mkchecksMod.F90
             mkdiagnosticsMod.F90
             mkdomainMod.F90
             mkesmfMod.F90
             mkfileMod.F90
             mkgdpMod.F90
             mkglacierregionMod.F90
             mkglcmecMod.F90
             mkharvestMod.F90
             mkindexmapMod.F90
             mkinputMod.F90
             mklaiMod.F90
             mklanwatMod.F90
             mkorganicMod.F90
             mkpeatMod.F90
             mkpioMod.F90
             mkpftMod.F90
             mkpftConstantsMod.F90
             mkpctPftTypeMod.F90
             mkpftUtilsMod.F90
             mksoilcolMod.F90
             mksoilfmaxMod.F90
             mksoiltexMod.F90
             mksoiltexnewMod.F90
             mksoildepthMod.F90
             mktopostatsMod.F90
             mkurbanparMod.F90
             mkutilsMod.F90
             mkvarctl.F90
             mkvarpar.F90
             mkvocefMod.F90
             mkVICparamsMod.F90
             nanMod.F90
             shr_const_mod.F90
             shr_kind_mod.F90
             shr_string_mod.F90
             shr_sys_mod.F90)

add_library(mksurfdata ${SRCFILES})

add_compile_definitions(PIO2)

message("CMAKE_CURRENT_BINARY_DIR is ${CMAKE_CURRENT_BINARY_DIR}")
message("PROJECT_BINARY_DIR is ${PROJECT_BINARY_DIR}")

target_include_directories (mksurfdata PRIVATE ${ESMF_F90COMPILEPATHS})
target_include_directories (mksurfdata PRIVATE /glade/work/jedwards/tools/pio/2.5.5/intel/19.1.1/mpt/2.22/include)
target_include_directories (mksurfdata PRIVATE /glade/u/apps/ch/opt/pnetcdf/1.12.2/mpt/2.22/intel/19.1.1//include)
target_include_directories (mksurfdata PRIVATE /glade/u/apps/ch/opt/netcdf-mpi/4.8.1/mpt/2.22/intel/19.1.1//include)

install(TARGETS mksurfdata)
