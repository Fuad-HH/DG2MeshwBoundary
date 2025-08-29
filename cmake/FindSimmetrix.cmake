# - Try to find SIMMETRIX libraries
# Once done this will define
#  SIMMETRIX_FOUND - System has SIMMETRIX
#  SIMMETRIX_INCLUDE_DIRS - The SIMMETRIX include directories
#  SIMMETRIX_LIBRARIES - The libraries needed to use SIMMETRIX
#  SIMMETRIX_DEFINITIONS - Compiler switches required for using SIMMETRIX
#
# This implementation assumes a SIMMETRIX install has the following structure
# VERSION/
#         include/*.h
#         lib/*.a

macro(simmetrixLibCheck libs isRequired)
    foreach(lib ${libs})
        unset(simmetrixlib CACHE)
        find_library(simmetrixlib "${lib}" PATHS ${SIMMETRIX_LIB_DIR})
        if(simmetrixlib MATCHES "^simmetrixlib-NOTFOUND$")
            if(${isRequired})
                message(FATAL_ERROR "SIMMETRIX library ${lib} not found in ${SIMMETRIX_LIB_DIR}")
            else()
                message("SIMMETRIX library ${lib} not found in ${SIMMETRIX_LIB_DIR}")
            endif()
        else()
            set("SIMMETRIX_${lib}_FOUND" TRUE CACHE INTERNAL "SIMMETRIX library present")
            set(SIMMETRIX_LIBS ${SIMMETRIX_LIBS} ${simmetrixlib})
        endif()
    endforeach()
endmacro(simmetrixLibCheck)

set(SIMMETRIX_LIBS "")

if (ENABLE_PPPL)
    set(SIMMETRIX_LIB_NAMES
            SimLicense  #-- valid for PPPL
            SimPartitionedMesh #-mpi
            SimMeshing
            SimMeshTools
            SimModel
            SimPartitionWrapper #-${SIM_MPI}
            SimAdvMeshing
            #SimField -- not valid for PPPL
    )
else()
    set(SIMMETRIX_LIB_NAMES
            SimDiscrete
            SimPartitionedMesh-mpi
            SimMeshing
            SimMeshTools
            SimModel
            SimPartitionWrapper-${SIM_MPI}
            SimAdvMeshing
            SimField
    )
endif()

simmetrixLibCheck("${SIMMETRIX_LIB_NAMES}" TRUE)


string(FIND "${SIMMETRIX_LIBS}" "/lib/" archStart)
string(FIND "${SIMMETRIX_LIBS}" "/libSim" archEnd)
math(EXPR archStart "${archStart}+5")
math(EXPR len "${archEnd}-${archStart}")
string(SUBSTRING "${SIMMETRIX_LIBS}" "${archStart}" "${len}" SIM_ARCHOS)
message(STATUS "SIM_ARCHOS ${SIM_ARCHOS}")

find_path(SIMMETRIX_INCLUDE_DIR
        NAMES MeshSim.h
        PATHS ${SIMMETRIX_INCLUDE_DIR})
if(NOT EXISTS "${SIMMETRIX_INCLUDE_DIR}")
    message(FATAL_ERROR "SIMMETRIX include dir not found")
endif()

string(REGEX REPLACE
        "/include$" ""
        SIMMETRIX_INSTALL_DIR
        "${SIMMETRIX_INCLUDE_DIR}")

set(SIMMETRIX_LIBRARIES ${SIMMETRIX_LIBS} )
set(SIMMETRIX_INCLUDE_DIRS ${SIMMETRIX_INCLUDE_DIR} )

if (SIM_ARCHOS STREQUAL x64_rhel8_gcc83)
    find_library(XDR_LIB tirpc)
    if(XDR_LIB)
        message(STATUS "Found XDR_LIB ${XDR_LIB}")
        set(SIMMETRIX_LIBS ${SIMMETRIX_LIBS} ${XDR_LIB})
    else()
        message(FATAL_ERROR "The libtirpc library was not found.  It defines xdr symbols "
                "(e.g., xdrmem_create) that are need by SimModSuite on systems using "
                "glibc newer than 2.32.  Note, glibc starting with 2.26 could optionally "
                "have been built without the xdr symbols.")
    endif()
endif()

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set PARMETIS_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(SIMMETRIX  DEFAULT_MSG
        SIMMETRIX_LIBS SIMMETRIX_INCLUDE_DIR)

mark_as_advanced(SIMMETRIX_INCLUDE_DIR SIMMETRIX_LIBS)

set(SIMMETRIX_LINK_LIBS "")
foreach(lib ${SIMMETRIX_LIB_NAMES})
    set(SIMMETRIX_LINK_LIBS "${SIMMETRIX_LINK_LIBS} -l${lib}")
endforeach()

#pkgconfig
set(prefix "${SIMMETRIX_INSTALL_DIR}")
set(includedir "${SIMMETRIX_INCLUDE_DIR}")
configure_file(
        "${CMAKE_HOME_DIRECTORY}/cmake/libSimmetrix.pc.in"
        "${CMAKE_BINARY_DIR}/libSimmetrix.pc"
        @ONLY)

INSTALL(FILES "${CMAKE_BINARY_DIR}/libSimmetrix.pc" DESTINATION lib/pkgconfig)