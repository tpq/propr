###############################################################################
# FindTestthat.cmake  -  Locate R's testthat package (headers & optional library)
# Usage:
#   find_package(Testthat [REQUIRED] [QUIET])
#
# It sets the following variables:
#   Testthat_FOUND           - True if testthat is found
#   Testthat_INCLUDE_DIRS    - Full path to testthat/include
#   Testthat_LIBRARY         - Path to libtestthat (if available)
#
# The following variables will be checked or set:
#   R_ROOT_PATH              - Optionally preset R_HOME path
###############################################################################

if (DEFINED ENV{R_HOME})
    set(R_ROOT_PATH "$ENV{R_HOME}")
else()
    execute_process(COMMAND R RHOME
            OUTPUT_VARIABLE _R_HOME
            ERROR_QUIET)
    string(REGEX REPLACE "\n" "" R_ROOT_PATH "${_R_HOME}")
endif()

execute_process(
        COMMAND ${R_ROOT_PATH}/bin/Rscript -e "cat(find.package('testthat'))"
        OUTPUT_VARIABLE Testthat_ROOT
        ERROR_QUIET
)
string(REGEX REPLACE "\n" "" Testthat_ROOT "${Testthat_ROOT}")

if (Testthat_ROOT)
    set(Testthat_INCLUDE_DIRS "${Testthat_ROOT}/include")
    find_path(
            _TESTTHAT_HDR
            NAMES "testthat.h"
            PATHS ${Testthat_INCLUDE_DIRS}
            NO_DEFAULT_PATH
    )
endif()

if (Testthat_ROOT)
    execute_process(
            COMMAND ${R_ROOT_PATH}/bin/Rscript -e "cat(find.package('testthat'))"
            OUTPUT_VARIABLE _TESTTHAT_LIB_PATH
            ERROR_QUIET
    )

    find_library(
            _TESTTHAT_LIB
            REQUIRED
            NAMES "testthat.so"
            PATHS ${_TESTTHAT_LIB_PATH}
            PATH_SUFFIXES "/libs" "/lib64" "/bin"
            NO_DEFAULT_PATH
    )
endif()

add_library(Testthat::Testthat UNKNOWN IMPORTED)

set(Testthat_VERSION_MAJOR 3)
set(Testthat_VERSION_MINOR 0)

include(FindPackageHandleStandardArgs)

set(Testthat_INCLUDE_DIRS ${_TESTTHAT_HDR})
set(Testthat_LIBRARY      ${_TESTTHAT_LIB})
find_package_handle_standard_args( Testthat DEFAULT_MSG Testthat_INCLUDE_DIRS Testthat_LIBRARY )

if (Testthat_FOUND)
    set_target_properties(Testthat::Testthat
            PROPERTIES
                INTERFACE_INCLUDE_DIRECTORIES
                    "${Testthat_INCLUDE_DIRS}"
                    IMPORTED_LOCATION  "${Testthat_LIBRARY}"
    )
    #if (Testthat_LIBRARY)
    #    set_target_properties(Testthat::Testthat
    #                PROPERTIES
    #                    INTERFACE_LINK_LIBRARIES
    #                        "${Testthat_LIBRARY}" )
    #endif()
    mark_as_advanced(Testthat_INCLUDE_DIRS Testthat_LIBRARY)
endif()
