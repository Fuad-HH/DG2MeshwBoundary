cmake_minimum_required(VERSION 3.0)
project(simmeshWboundaryNodes)

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(SimModSuite REQUIRED)
include_directories(${SimModSuite_INCLUDE_DIRS})
find_package(Simmetrix MODULE REQUIRED)
include_directories(${SIMMETRIX_INCLUDE_DIR})

add_executable(simmeshWboundaryNodes simmeshWboundaryNodes.cpp)
target_link_libraries(simmeshWboundaryNodes PRIVATE ${SIMMETRIX_LIBS} ${SIMMODSUITE_LIBS})
target_include_directories(simmeshWboundaryNodes PUBLIC ${SIMMETRIX_INCLUDE_DIR})

set_target_properties(simmeshWboundaryNodes PROPERTIES
        INSTALL_RPATH_USE_LINK_PATH TRUE
        INSTALL_RPATH "$ORIGIN/../lib:$ORIGIN/../lib64"
)

install(TARGETS simmeshWboundaryNodes
         RUNTIME DESTINATION bin)
