# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file LICENSE.rst or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION ${CMAKE_VERSION}) # this file comes with cmake

# If CMAKE_DISABLE_SOURCE_CHANGES is set to true and the source directory is an
# existing directory in our source tree, calling file(MAKE_DIRECTORY) on it
# would cause a fatal error, even though it would be a no-op.
if(NOT EXISTS "/Users/michel/Desktop/Geomatics/Year_2/Thesis/C++/external/masbcpp")
  file(MAKE_DIRECTORY "/Users/michel/Desktop/Geomatics/Year_2/Thesis/C++/external/masbcpp")
endif()
file(MAKE_DIRECTORY
  "/Users/michel/Desktop/Geomatics/Year_2/Thesis/C++/build-test/_deps/masbcpp-build"
  "/Users/michel/Desktop/Geomatics/Year_2/Thesis/C++/build-test/masbcpp_ext-prefix"
  "/Users/michel/Desktop/Geomatics/Year_2/Thesis/C++/build-test/masbcpp_ext-prefix/tmp"
  "/Users/michel/Desktop/Geomatics/Year_2/Thesis/C++/build-test/masbcpp_ext-prefix/src/masbcpp_ext-stamp"
  "/Users/michel/Desktop/Geomatics/Year_2/Thesis/C++/build-test/masbcpp_ext-prefix/src"
  "/Users/michel/Desktop/Geomatics/Year_2/Thesis/C++/build-test/masbcpp_ext-prefix/src/masbcpp_ext-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "/Users/michel/Desktop/Geomatics/Year_2/Thesis/C++/build-test/masbcpp_ext-prefix/src/masbcpp_ext-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "/Users/michel/Desktop/Geomatics/Year_2/Thesis/C++/build-test/masbcpp_ext-prefix/src/masbcpp_ext-stamp${cfgdir}") # cfgdir has leading slash
endif()
