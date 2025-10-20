cmake_minimum_required(VERSION 3.0)
project(DG2MeshwBoundary)

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(SimModSuite REQUIRED)
include_directories(${SimModSuite_INCLUDE_DIR})

add_executable(dg2MeshwBoundary DG2MeshwBoundary.cpp MeshConfig.cpp)
target_link_libraries(dg2MeshwBoundary ${SIMMODSUITE_LIBS})
target_include_directories(dg2MeshwBoundary PRIVATE ${SIMMODSUITE_INCLUDE_DIR} ${CMAKE_CURRENT_SOURCE_DIR})

set_target_properties(dg2MeshwBoundary PROPERTIES
        INSTALL_RPATH_USE_LINK_PATH TRUE
        INSTALL_RPATH "$ORIGIN/../lib:$ORIGIN/../lib64"
)

install(TARGETS dg2MeshwBoundary
         RUNTIME DESTINATION bin)
