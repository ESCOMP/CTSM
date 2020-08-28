# 
# Author: Ali Samii - The University of Texas at Austin
#
# Distributed under GPL2. For more info refer to:
#    https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html
#
# 
# FindESMF
# --------
# 
# This script tries to find the ESMF library. You have to define
# the path to esmf.mk file in your installation directory.
#
# There are plans to extend this script to find ESMF automatically,
# but until then, you should set the environment variable
#
#      ESMF_CONFIG_FILE = /path/to/esmf.mk
#
# in your installation directory. The output will be
#
#      ESMF_LINK_LINE       : All the libraries and link line stuff 
#      ESMF_COMPILER_LINE   : All the compiler flags and include dirs
#
#

# Defining the ${Esc} for syntax coloring.
string(ASCII 27 Esc)

message ("Parsing ESMF_CONFIG_FILE: " $ENV{ESMF_CONFIG_FILE})

file(STRINGS "$ENV{ESMF_CONFIG_FILE}" all_vars)
foreach(str ${all_vars})
  string(REGEX MATCH "^[^#]" def ${str})
  if (def)
    string(REGEX MATCH "^[^=]+" var_name ${str})
    string(REGEX MATCH "=(.+)$" var_def ${str})
    set(var_def ${CMAKE_MATCH_1})
    set(${var_name} ${var_def})
    mark_as_advanced (${var_name})
  endif()
endforeach()

set (ESMF_LINK_LINE "${ESMF_F90LINKOPTS} \
                     ${ESMF_F90LINKRPATHS} \
                     ${ESMF_F90LINKPATHS} \
                     ${ESMF_F90ESMFLINKLIBS}")

set (ESMF_COMPILER_LINE "${ESMF_F90COMPILEOPTS} \
                         ${ESMF_F90COMPILEPATHS} \
                         ${ESMF_F90COMPILEFREENOCPP} \
                         ${ESMF_CXXCOMPILEPATHS}")

