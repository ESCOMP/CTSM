# (C) Copyright 2011- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

# Try to find NetCDF includes and library.
# Supports static and shared libaries and allows each component to be found in sepearte prefixes.
#
# This module defines
#
#   - NetCDF_FOUND                - System has NetCDF
#   - NetCDF_INCLUDE_DIRS         - the NetCDF include directories
#   - NetCDF_VERSION              - the version of NetCDF
#   - NetCDF_CONFIG_EXECUTABLE    - the netcdf-config executable if found
#   - NetCDF_PARALLEL             - Boolean True if NetCDF4 has parallel IO support via hdf5 and/or pnetcdf
#   - NetCDF_HAS_PNETCDF          - Boolean True if NetCDF4 has pnetcdf support
#
# Deprecated Defines
#   - NetCDF_LIBRARIES            - [Deprecated] Use NetCDF::NetCDF_<LANG> targets instead.
#
#
# Following components are available:
#
#   - C                           - C interface to NetCDF          (netcdf)
#   - CXX                         - CXX4 interface to NetCDF       (netcdf_c++4)
#   - Fortran                     - Fortran interface to NetCDF    (netcdff)
#
# For each component the following are defined:
#
#   - NetCDF_<comp>_FOUND         - whether the component is found
#   - NetCDF_<comp>_LIBRARIES     - the libraries for the component
#   - NetCDF_<comp>_LIBRARY_SHARED - Boolean is true if libraries for component are shared
#   - NetCDF_<comp>_INCLUDE_DIRS  - the include directories for specified component
#   - NetCDF::NetCDF_<comp>       - target of component to be used with target_link_libraries()
#
# The following paths will be searched in order if set in CMake (first priority) or environment (second priority)
#
#   - NetCDF_ROOT                 - root of NetCDF installation
#   - NetCDF_PATH                 - root of NetCDF installation
#
# The search process begins with locating NetCDF Include headers.  If these are in a non-standard location,
# set one of the following CMake or environment variables to point to the location:
#
#  - NetCDF_INCLUDE_DIR or NetCDF_${comp}_INCLUDE_DIR
#  - NetCDF_INCLUDE_DIRS or NetCDF_${comp}_INCLUDE_DIR
#
# Notes:
#
#   - Use "NetCDF::NetCDF_<LANG>" targets only.  NetCDF_LIBRARIES exists for backwards compatibility and should not be used.
#     - These targets have all the knowledge of include directories and library search directories, and a single
#       call to target_link_libraries will provide all these transitive properties to your target.  Normally all that is
#       needed to build and link against NetCDF is, e.g.:
#           target_link_libraries(my_c_tgt PUBLIC NetCDF::NetCDF_C)
#   - "NetCDF" is always the preferred naming for this package, its targets, variables, and environment variables
#     - For compatibility, some variables are also set/checked using alternate names NetCDF4, NETCDF, or NETCDF4
#     - Environments relying on these older environment variable names should move to using a "NetCDF_ROOT" environment variable
#   - Preferred component capitalization follows the CMake LANGUAGES variables: i.e., C, Fortran, CXX
#     - For compatibility, alternate capitalizations are supported but should not be used.
#   - If no components are defined, all components will be searched
#

list( APPEND _possible_components C CXX Fortran )

## Include names for each component
set( NetCDF_C_INCLUDE_NAME          netcdf.h )
set( NetCDF_CXX_INCLUDE_NAME        netcdf )
set( NetCDF_Fortran_INCLUDE_NAME    netcdf.mod )

## Library names for each component
set( NetCDF_C_LIBRARY_NAME          netcdf )
set( NetCDF_CXX_LIBRARY_NAME        netcdf_c++4 )
set( NetCDF_Fortran_LIBRARY_NAME    netcdff )

## Enumerate search components
foreach( _comp ${_possible_components} )
  string( TOUPPER "${_comp}" _COMP )
  set( _arg_${_COMP} ${_comp} )
  set( _name_${_COMP} ${_comp} )
endforeach()

set( _search_components C)
foreach( _comp ${${CMAKE_FIND_PACKAGE_NAME}_FIND_COMPONENTS} )
  string( TOUPPER "${_comp}" _COMP )
  set( _arg_${_COMP} ${_comp} )
  list( APPEND _search_components ${_name_${_COMP}} )
  if( NOT _name_${_COMP} )
    message(SEND_ERROR "Find${CMAKE_FIND_PACKAGE_NAME}: COMPONENT ${_comp} is not a valid component. Valid components: ${_possible_components}" )
  endif()
endforeach()
list( REMOVE_DUPLICATES _search_components )

## Search hints for finding include directories and libraries
foreach( _comp IN ITEMS "_" "_C_" "_Fortran_" "_CXX_" )
  foreach( _name IN ITEMS NetCDF4 NetCDF NETCDF4 NETCDF )
    foreach( _var IN ITEMS ROOT PATH )
      list(APPEND _search_hints ${${_name}${_comp}${_var}} $ENV{${_name}${_comp}${_var}} )
      list(APPEND _include_search_hints
                ${${_name}${_comp}INCLUDE_DIR} $ENV{${_name}${_comp}INCLUDE_DIR}
                ${${_name}${_comp}INCLUDE_DIRS} $ENV{${_name}${_comp}INCLUDE_DIRS} )
    endforeach()
  endforeach()
endforeach()
#Old-school HPC module env variable names
foreach( _name IN ITEMS NetCDF4 NetCDF NETCDF4 NETCDF )
  foreach( _comp IN ITEMS "_C" "_Fortran" "_CXX" )
    list(APPEND _search_hints ${${_name}}         $ENV{${_name}})
    list(APPEND _search_hints ${${_name}${_comp}} $ENV{${_name}${_comp}})
  endforeach()
endforeach()

## Find headers for each component
set(NetCDF_INCLUDE_DIRS)
set(_new_search_components)
foreach( _comp IN LISTS _search_components )
  if(NOT ${PROJECT_NAME}_NetCDF_${_comp}_FOUND)
      list(APPEND _new_search_components ${_comp})
  endif()
  find_file(NetCDF_${_comp}_INCLUDE_FILE
    NAMES ${NetCDF_${_comp}_INCLUDE_NAME}
    DOC "NetCDF ${_comp} include directory"
    HINTS ${_include_search_hints} ${_search_hints}
    PATH_SUFFIXES include include/netcdf
  )
  mark_as_advanced(NetCDF_${_comp}_INCLUDE_FILE)
  message(DEBUG "NetCDF_${_comp}_INCLUDE_FILE: ${NetCDF_${_comp}_INCLUDE_FILE}")
  if( NetCDF_${_comp}_INCLUDE_FILE )
    get_filename_component(NetCDF_${_comp}_INCLUDE_FILE ${NetCDF_${_comp}_INCLUDE_FILE} ABSOLUTE)
    get_filename_component(NetCDF_${_comp}_INCLUDE_DIR ${NetCDF_${_comp}_INCLUDE_FILE} DIRECTORY)
    list(APPEND NetCDF_INCLUDE_DIRS ${NetCDF_${_comp}_INCLUDE_DIR})
  endif()
endforeach()
if(NetCDF_INCLUDE_DIRS)
    list(REMOVE_DUPLICATES NetCDF_INCLUDE_DIRS)
endif()
set(NetCDF_INCLUDE_DIRS "${NetCDF_INCLUDE_DIRS}" CACHE STRING "NetCDF Include directory paths" FORCE)

## Find n*-config executables for search components
foreach( _comp IN LISTS _search_components )
  if( _comp MATCHES "^(C)$" )
    set(_conf "c")
  elseif( _comp MATCHES "^(Fortran)$" )
    set(_conf "f")
  elseif( _comp MATCHES "^(CXX)$" )
    set(_conf "cxx4")
  endif()
  find_program( NetCDF_${_comp}_CONFIG_EXECUTABLE
      NAMES n${_conf}-config
    HINTS ${NetCDF_INCLUDE_DIRS} ${_include_search_hints} ${_search_hints}
    PATH_SUFFIXES bin Bin ../bin ../../bin
      DOC "NetCDF n${_conf}-config helper" )
    message(DEBUG "NetCDF_${_comp}_CONFIG_EXECUTABLE: ${NetCDF_${_comp}_CONFIG_EXECUTABLE}")
endforeach()

set(_C_libs_flag --libs)
set(_Fortran_libs_flag --flibs)
set(_CXX_libs_flag --libs)
set(_C_includes_flag --includedir)
set(_Fortran_includes_flag --includedir)
set(_CXX_includes_flag --includedir)
function(netcdf_config exec flag output_var)
  set(${output_var} False PARENT_SCOPE)
  if( exec )
    execute_process( COMMAND ${exec} ${flag} RESULT_VARIABLE _ret OUTPUT_VARIABLE _val)
    if( _ret EQUAL 0 )
      string( STRIP ${_val} _val )
      set( ${output_var} ${_val} PARENT_SCOPE )
    endif()
  endif()
endfunction()

## Detect additional package properties
netcdf_config(${NetCDF_C_CONFIG_EXECUTABLE} --has-parallel4 _val)
if( NOT _val MATCHES "^(yes|no)$" )
  netcdf_config(${NetCDF_C_CONFIG_EXECUTABLE} --has-parallel _val)
endif()
if( _val MATCHES "^(yes)$" )
  set(NetCDF_PARALLEL TRUE CACHE STRING "NetCDF has parallel IO capability via pnetcdf or hdf5." FORCE)
else()
  set(NetCDF_PARALLEL FALSE CACHE STRING "NetCDF has no parallel IO capability." FORCE)
endif()

if(NetCDF_PARALLEL)
  find_package(MPI REQUIRED)
endif()

## Find libraries for each component
set( NetCDF_LIBRARIES )
foreach( _comp IN LISTS _search_components )
  string( TOUPPER "${_comp}" _COMP )

  find_library( NetCDF_${_comp}_LIBRARY
    NAMES ${NetCDF_${_comp}_LIBRARY_NAME}
    DOC "NetCDF ${_comp} library"
    HINTS ${NetCDF_${_comp}_INCLUDE_DIRS} ${_search_hints}
    PATH_SUFFIXES lib64 lib ../lib64 ../lib ../../lib64 ../../lib )
  mark_as_advanced( NetCDF_${_comp}_LIBRARY )
  get_filename_component(NetCDF_${_comp}_LIBRARY ${NetCDF_${_comp}_LIBRARY} ABSOLUTE)
  set(NetCDF_${_comp}_LIBRARY ${NetCDF_${_comp}_LIBRARY} CACHE STRING "NetCDF ${_comp} library" FORCE)
  message(DEBUG "NetCDF_${_comp}_LIBRARY: ${NetCDF_${_comp}_LIBRARY}")

  if( NetCDF_${_comp}_LIBRARY )
    if( NetCDF_${_comp}_LIBRARY MATCHES ".a$" )
      set( NetCDF_${_comp}_LIBRARY_SHARED FALSE )
      set( _library_type STATIC)
    else()
      list( APPEND NetCDF_LIBRARIES ${NetCDF_${_comp}_LIBRARY} )
      set( NetCDF_${_comp}_LIBRARY_SHARED TRUE )
      set( _library_type SHARED)
    endif()
  endif()

  #Use nc-config to set per-component LIBRARIES variable if possible
  netcdf_config( ${NetCDF_${_comp}_CONFIG_EXECUTABLE} ${_${_comp}_libs_flag} _val )
  if( _val )
    set( NetCDF_${_comp}_LIBRARIES ${_val} )
    if(NOT NetCDF_${_comp}_LIBRARY_SHARED AND NOT NetCDF_${_comp}_FOUND) #Static targets should use nc_config to get a proper link line with all necessary static targets.
      list( APPEND NetCDF_LIBRARIES ${NetCDF_${_comp}_LIBRARIES} )
    endif()
  else()
    set( NetCDF_${_comp}_LIBRARIES ${NetCDF_${_comp}_LIBRARY} )
    if(NOT NetCDF_${_comp}_LIBRARY_SHARED)
      message(SEND_ERROR "Unable to properly find NetCDF.  Found static libraries at: ${NetCDF_${_comp}_LIBRARY} but could not run nc-config: ${NetCDF_CONFIG_EXECUTABLE}")
    endif()
  endif()

  #Use nc-config to set per-component INCLUDE_DIRS variable if possible
  netcdf_config( ${NetCDF_${_comp}_CONFIG_EXECUTABLE} ${_${_comp}_includes_flag} _val )
  if( _val )
    string( REPLACE " " ";" _val ${_val} )
    set( NetCDF_${_comp}_INCLUDE_DIRS ${_val} )
  else()
    set( NetCDF_${_comp}_INCLUDE_DIRS ${NetCDF_${_comp}_INCLUDE_DIR} )
  endif()

  if( NetCDF_${_comp}_LIBRARIES AND NetCDF_${_comp}_INCLUDE_DIRS )
    set( ${CMAKE_FIND_PACKAGE_NAME}_${_arg_${_COMP}}_FOUND TRUE )
    if (NOT TARGET NetCDF::NetCDF_${_comp})
      add_library(NetCDF::NetCDF_${_comp} ${_library_type} IMPORTED)
      set_target_properties(NetCDF::NetCDF_${_comp} PROPERTIES
        IMPORTED_LOCATION ${NetCDF_${_comp}_LIBRARY}
        INTERFACE_INCLUDE_DIRECTORIES "${NetCDF_${_comp}_INCLUDE_DIRS}"
        INTERFACE_LINK_LIBRARIES ${NetCDF_${_comp}_LIBRARIES} )
      if( NOT _comp MATCHES "^(C)$" )
        target_link_libraries(NetCDF::NetCDF_${_comp} INTERFACE NetCDF::NetCDF_C)
      endif()
      if(MPI_${_comp}_FOUND)
        target_link_libraries(NetCDF::NetCDF_${_comp} INTERFACE MPI::MPI_${_comp})
      endif()
    endif()
  endif()
endforeach()
if(NetCDF_LIBRARIES AND NetCDF_${_comp}_LIBRARY_SHARED)
    list(REMOVE_DUPLICATES NetCDF_LIBRARIES)
endif()
set(NetCDF_LIBRARIES "${NetCDF_LIBRARIES}" CACHE STRING "NetCDF library targets" FORCE)

## Find version via netcdf-config if possible
if (NetCDF_INCLUDE_DIRS)
  if( NetCDF_C_CONFIG_EXECUTABLE )
    netcdf_config( ${NetCDF_C_CONFIG_EXECUTABLE} --version _vers )
    if( _vers )
      string(REGEX REPLACE ".* ((([0-9]+)\\.)+([0-9]+)).*" "\\1" NetCDF_VERSION "${_vers}" )
    endif()
  else()
    foreach( _dir IN LISTS NetCDF_INCLUDE_DIRS)
      if( EXISTS "${_dir}/netcdf_meta.h" )
        file(STRINGS "${_dir}/netcdf_meta.h" _netcdf_version_lines
        REGEX "#define[ \t]+NC_VERSION_(MAJOR|MINOR|PATCH|NOTE)")
        string(REGEX REPLACE ".*NC_VERSION_MAJOR *\([0-9]*\).*" "\\1" _netcdf_version_major "${_netcdf_version_lines}")
        string(REGEX REPLACE ".*NC_VERSION_MINOR *\([0-9]*\).*" "\\1" _netcdf_version_minor "${_netcdf_version_lines}")
        string(REGEX REPLACE ".*NC_VERSION_PATCH *\([0-9]*\).*" "\\1" _netcdf_version_patch "${_netcdf_version_lines}")
        string(REGEX REPLACE ".*NC_VERSION_NOTE *\"\([^\"]*\)\".*" "\\1" _netcdf_version_note "${_netcdf_version_lines}")
        set(NetCDF_VERSION "${_netcdf_version_major}.${_netcdf_version_minor}.${_netcdf_version_patch}${_netcdf_version_note}")
        unset(_netcdf_version_major)
        unset(_netcdf_version_minor)
        unset(_netcdf_version_patch)
        unset(_netcdf_version_note)
        unset(_netcdf_version_lines)
      endif()
    endforeach()
  endif()
endif ()

## Finalize find_package
include(FindPackageHandleStandardArgs)

if(NOT NetCDF_FOUND OR _new_search_components)
    find_package_handle_standard_args( ${CMAKE_FIND_PACKAGE_NAME}
        REQUIRED_VARS NetCDF_INCLUDE_DIRS NetCDF_LIBRARIES
        VERSION_VAR NetCDF_VERSION
        HANDLE_COMPONENTS )
endif()

foreach( _comp IN LISTS _search_components )
    if( NetCDF_${_comp}_FOUND )
        #Record found components to avoid duplication in NetCDF_LIBRARIES for static libraries
        set(NetCDF_${_comp}_FOUND ${NetCDF_${_comp}_FOUND} CACHE BOOL "NetCDF ${_comp} Found" FORCE)
        #Set a per-package, per-component found variable to communicate between multiple calls to find_package()
        set(${PROJECT_NAME}_NetCDF_${_comp}_FOUND True)
    endif()
endforeach()

if( ${CMAKE_FIND_PACKAGE_NAME}_FOUND AND NOT ${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY AND _new_search_components)
  message( STATUS "Find${CMAKE_FIND_PACKAGE_NAME} defines targets:" )
  message( STATUS "  - NetCDF_VERSION [${NetCDF_VERSION}]")
  message( STATUS "  - NetCDF_PARALLEL [${NetCDF_PARALLEL}]")
  foreach( _comp IN LISTS _new_search_components )
    string( TOUPPER "${_comp}" _COMP )
    message( STATUS "  - NetCDF_${_comp}_CONFIG_EXECUTABLE [${NetCDF_${_comp}_CONFIG_EXECUTABLE}]")
    if( ${CMAKE_FIND_PACKAGE_NAME}_${_arg_${_COMP}}_FOUND )
      get_filename_component(_root ${NetCDF_${_comp}_INCLUDE_DIR}/.. ABSOLUTE)
      if( NetCDF_${_comp}_LIBRARY_SHARED )
        message( STATUS "  - NetCDF::NetCDF_${_comp} [SHARED] [Root: ${_root}] Lib: ${NetCDF_${_comp}_LIBRARY} ")
      else()
        message( STATUS "  - NetCDF::NetCDF_${_comp} [STATIC] [Root: ${_root}] Lib: ${NetCDF_${_comp}_LIBRARY} ")
      endif()
    endif()
  endforeach()
endif()

foreach( _prefix NetCDF NetCDF4 NETCDF NETCDF4 ${CMAKE_FIND_PACKAGE_NAME} )
  set( ${_prefix}_INCLUDE_DIRS ${NetCDF_INCLUDE_DIRS} )
  set( ${_prefix}_LIBRARIES    ${NetCDF_LIBRARIES})
  set( ${_prefix}_VERSION      ${NetCDF_VERSION} )
  set( ${_prefix}_FOUND        ${${CMAKE_FIND_PACKAGE_NAME}_FOUND} )
  set( ${_prefix}_CONFIG_EXECUTABLE ${NetCDF_CONFIG_EXECUTABLE} )
  set( ${_prefix}_PARALLEL ${NetCDF_PARALLEL} )

  foreach( _comp ${_search_components} )
    string( TOUPPER "${_comp}" _COMP )
    set( _arg_comp ${_arg_${_COMP}} )
    set( ${_prefix}_${_comp}_FOUND     ${${CMAKE_FIND_PACKAGE_NAME}_${_arg_comp}_FOUND} )
    set( ${_prefix}_${_COMP}_FOUND     ${${CMAKE_FIND_PACKAGE_NAME}_${_arg_comp}_FOUND} )
    set( ${_prefix}_${_arg_comp}_FOUND ${${CMAKE_FIND_PACKAGE_NAME}_${_arg_comp}_FOUND} )

    set( ${_prefix}_${_comp}_LIBRARIES     ${NetCDF_${_comp}_LIBRARIES} )
    set( ${_prefix}_${_COMP}_LIBRARIES     ${NetCDF_${_comp}_LIBRARIES} )
    set( ${_prefix}_${_arg_comp}_LIBRARIES ${NetCDF_${_comp}_LIBRARIES} )

    set( ${_prefix}_${_comp}_INCLUDE_DIRS     ${NetCDF_${_comp}_INCLUDE_DIRS} )
    set( ${_prefix}_${_COMP}_INCLUDE_DIRS     ${NetCDF_${_comp}_INCLUDE_DIRS} )
    set( ${_prefix}_${_arg_comp}_INCLUDE_DIRS ${NetCDF_${_comp}_INCLUDE_DIRS} )
  endforeach()
endforeach()