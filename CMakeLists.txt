cmake_minimum_required(VERSION 3.16)
project(Thesis)

set(CMAKE_CXX_STANDARD 20)

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


target_link_libraries(alpha_wrap CGAL::CGAL)

#val3dity
find_program(VAL3DITY_EXECUTABLE val3dity /opt/homebrew/bin)

if (NOT VAL3DITY_EXECUTABLE)
    message(FATAL_ERROR "val3dity not found at /opt/homebrew/bin")
endif()

add_definitions(-DVAL3DITY_PATH="${VAL3DITY_EXECUTABLE}")
