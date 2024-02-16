cmake_minimum_required(VERSION 3.10)
project(mksurfdata Fortran)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../cmake")
find_package(NetCDF 4.7.4 REQUIRED Fortran)
find_package(ESMF 8.2.0 REQUIRED )

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
             mkpeatMod.F90
             mkpioMod.F90
             mkpftMod.F90
             mkpftConstantsMod.F90
             mkpctPftTypeMod.F90
             mkpftUtilsMod.F90
             mksoilcolMod.F90
             mksoilfmaxMod.F90
             mksoiltexMod.F90
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
             shr_sys_mod.F90
             mksurfdata.F90)

add_compile_definitions(PIO2)

add_library(pnetcdf STATIC IMPORTED)
set_property(TARGET pnetcdf PROPERTY IMPORTED_LOCATION $ENV{PNETCDF})

add_library(pioc STATIC IMPORTED)
add_library(piof STATIC IMPORTED)
set_property(TARGET pioc PROPERTY IMPORTED_LOCATION $ENV{PIO}/lib/libpioc.so)
set_property(TARGET piof PROPERTY IMPORTED_LOCATION $ENV{PIO}/lib/libpiof.so)
add_executable(mksurfdata ${SRCFILES})
target_link_libraries(mksurfdata PRIVATE esmf piof pioc)
target_include_directories (mksurfdata PRIVATE ${ESMF_F90COMPILEPATHS})
target_include_directories (mksurfdata PRIVATE $ENV{PIO}/include)
target_include_directories (mksurfdata PRIVATE ${PNETCDF}/include)
target_include_directories (mksurfdata PRIVATE ${NETCDF}/include)

install(TARGETS mksurfdata)
