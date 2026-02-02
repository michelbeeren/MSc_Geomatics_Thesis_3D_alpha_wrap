cmake_minimum_required(VERSION 3.16)
project(Thesis)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Use local modified CGAL headers first, can exclude this, than normal alpha wrap algorithm is used
#include_directories(
#        ${CMAKE_SOURCE_DIR}/cgal_override
#)

find_package(CGAL REQUIRED)

add_executable(alpha_wrap src/main.cpp
        src/val3dity.cpp
        src/val3dity.h
        src/alpha_wrap.cpp
        src/alpha_wrap.h
        src/edge_refinement.cpp
        src/edge_refinement.h
        src/octree.cpp
        src/octree.h
        src/MAT.cpp
        src/MAT.h
)


target_link_libraries(alpha_wrap PRIVATE CGAL::CGAL)

#val3dity
find_program(VAL3DITY_EXECUTABLE val3dity /opt/homebrew/bin)

if (NOT VAL3DITY_EXECUTABLE)
    message(FATAL_ERROR "val3dity not found at /opt/homebrew/bin")
endif()

target_compile_definitions(alpha_wrap PRIVATE VAL3DITY_PATH="${VAL3DITY_EXECUTABLE}")

# ---------------- masbcpp (ExternalProject) ----------------
include(ExternalProject)

set(MASBCPP_SOURCE_DIR ${CMAKE_SOURCE_DIR}/external/masbcpp)
set(MASBCPP_BINARY_DIR ${CMAKE_BINARY_DIR}/_deps/masbcpp-build)

ExternalProject_Add(masbcpp_ext
        SOURCE_DIR ${MASBCPP_SOURCE_DIR}
        BINARY_DIR ${MASBCPP_BINARY_DIR}
        CMAKE_ARGS
        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
        -DCMAKE_PREFIX_PATH=/opt/homebrew
        BUILD_COMMAND ${CMAKE_COMMAND} --build .
        INSTALL_COMMAND ""
)

# Build masbcpp whenever alpha_wrap is built
add_dependencies(alpha_wrap masbcpp_ext)

# Provide the paths to the built executables to your C++ code
target_compile_definitions(alpha_wrap PRIVATE
        MASBCPP_COMPUTE_MA_PATH="${MASBCPP_BINARY_DIR}/compute_ma"
        MASBCPP_SIMPLIFY_PATH="${MASBCPP_BINARY_DIR}/simplify"
)

target_include_directories(alpha_wrap PRIVATE
        ${CMAKE_SOURCE_DIR}/external/masbcpp/thirdparty/cnpy
)

target_sources(alpha_wrap PRIVATE
        ${CMAKE_SOURCE_DIR}/external/masbcpp/thirdparty/cnpy/cnpy.cpp
)
