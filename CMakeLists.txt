cmake_minimum_required(VERSION 3.16)
project(Thesis)

set(CMAKE_CXX_STANDARD 20)

find_package(CGAL REQUIRED)

add_executable(alpha_wrap_example 3D_alpha_wrap.cpp)
add_executable(main_example main.cpp)


target_link_libraries(alpha_wrap_example CGAL::CGAL)

#val3dity
find_program(VAL3DITY_EXECUTABLE val3dity /opt/homebrew/bin)

if (NOT VAL3DITY_EXECUTABLE)
    message(FATAL_ERROR "val3dity not found at /opt/homebrew/bin")
endif()

add_definitions(-DVAL3DITY_PATH="${VAL3DITY_EXECUTABLE}")
