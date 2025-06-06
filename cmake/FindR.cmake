##########################################################################
# Copyright (c) 2023, King Abdullah University of Science and Technology
# All rights reserved.
# MPCR is an R package provided by the STSDS group at KAUST
##########################################################################

# - Find the R and Rcpp library
#
# Usage:
#   find_package(R [REQUIRED] [QUIET] )
#
# It sets the following variables:
#   R _FOUND               ... true if R and Rcpp is found on the system
#   R_LIB                  ... full path to R library
#   R_INCLUDE              ... R and Rcpp include directory
#
# The following variables will be checked by the function
#   R_ROOT_PATH          ... if set, the R libraries are exclusively searched
#                               under this path
#   R_LIB_PATH          ... if set, the Rcpp libraries are exclusively searched
#                              under this path
#   R_INCLUDE_DIR       ... R headers path from the R-environment variables.

if (DEFINED ENV{R_HOME})
    set(R_ROOT_PATH "$ENV{R_HOME}")

else ()
    execute_process(COMMAND R RHOME OUTPUT_VARIABLE R_HOME)
    string(REGEX REPLACE "\n" "" R_HOME "${R_HOME}")
    set(R_ROOT_PATH "${R_HOME}")

endif ()

if (NOT R_INCLUDE_PATH)
    execute_process(COMMAND ${R_ROOT_PATH}/bin/Rscript -e "cat(Sys.getenv('R_INCLUDE_DIR'))" OUTPUT_VARIABLE R_INCLUDE_DIR)
    string(REGEX REPLACE "\n" "" R_INCLUDE_DIR "${R_INCLUDE_DIR}")
    set(R_INCLUDE_PATH "${R_INCLUDE_DIR}")
    message("R Include Path :  " ${R_INCLUDE_PATH})
endif ()


if (NOT RCPP_LIB_PATH)
    execute_process(COMMAND ${R_ROOT_PATH}/bin/Rscript ${CMAKE_MODULE_PATH}/FindRLibraryPath.R OUTPUT_VARIABLE RCPP_LIB_PATH)
    set(RCPP_LIB_PATH ${RCPP_LIB_PATH})
    message("Rcpp Lib Path :  " ${RCPP_LIB_PATH})
endif ()

message("R Home Path :  " ${R_ROOT_PATH})

if (R_ROOT_PATH)

    if (APPLE)
        find_library(
                R_DYN_LIB
                NAMES "libR.dylib"
                PATHS ${R_ROOT_PATH}
                PATH_SUFFIXES "lib" "lib64" "bin"
                NO_DEFAULT_PATH
        )

    else ()
        #find libs
        find_library(
                R_DYN_LIB
                NAMES "libR.so"
                PATHS ${R_ROOT_PATH}
                PATH_SUFFIXES "lib" "lib64" "bin"
                NO_DEFAULT_PATH
        )

    endif ()

    if (R_DYN_LIB MATCHES R_DYN_LIB-NOTFOUND)
        set(R_DYN_LIB "")
        message("R is built with no dynamic library support")
    endif ()

    set(R_LIB
            ${R_DYN_LIB}
    )

else ()
    error("R is not installed ")
endif (R_ROOT_PATH)


if (R_INCLUDE_PATH)
    # find includes
    find_path(
            R_INCLUDE_DIRS
            REQUIRED
            NAMES "R.h"
            PATHS ${R_INCLUDE_PATH}
            PATH_SUFFIXES "include"
            NO_DEFAULT_PATH
    )
endif ()

if (RCPP_LIB_PATH)
    #find libs
    find_library(
            RCPP_LIB
            REQUIRED
            NAMES "Rcpp.so"
            PATHS ${RCPP_LIB_PATH}
            PATH_SUFFIXES "/libs" "/lib64" "/bin"
            NO_DEFAULT_PATH
    )

    # find includes
    find_path(
            RCPP_INCLUDE_DIRS
            REQUIRED
            NAMES "Rcpp.h"
            PATHS ${RCPP_LIB_PATH}
            PATH_SUFFIXES "/include"
            NO_DEFAULT_PATH
    )


else ()
    message("Rcpp is not installed ...")
endif (RCPP_LIB_PATH)

set(R_INCLUDE
        ${R_INCLUDE}
        ${R_INCLUDE_DIRS}
        ${RCPP_INCLUDE_DIRS}
)

add_library(R INTERFACE IMPORTED)

if (R_LIB)
    include(FindPackageHandleStandardArgs)
    find_package_handle_standard_args(R DEFAULT_MSG
            R_INCLUDE R_LIB)

    include_directories(${R_INCLUDE})
    mark_as_advanced(R_INCLUDE R_INCLUDE_DIRS RCPP_INCLUDE_DIRS R_LIB RCPP_LIB)

else ()
    set_target_properties(R
            PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${R_INCLUDE}"
            IMPORTED_LOCATION ${RCPP_LIB}
    )
    include(FindPackageHandleStandardArgs)
    find_package_handle_standard_args(R DEFAULT_MSG
            R_INCLUDE)

    include_directories(${R_INCLUDE})
    mark_as_advanced(R_INCLUDE R_INCLUDE_DIRS RCPP_INCLUDE_DIRS RCPP_LIB)
endif ()