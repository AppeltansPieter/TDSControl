cmake_minimum_required(VERSION 3.16.0)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
option(ENABLE_WARNINGS_AS_ERRORS "Treat warnings as errors" OFF)

file(READ "${CMAKE_CURRENT_SOURCE_DIR}/VERSION" VERSION_CONTENTS)
string(REGEX MATCH "VERSION = \"([0-9]+\.[0-9]+\.[0-9])\"" _
             ${VERSION_CONTENTS})
set(PROJECT_VERSION "${CMAKE_MATCH_1}")

project(
  tds-control
  VERSION ${PROJECT_VERSION}
  LANGUAGES C CXX)

include(FetchContent)

set(TDS_CONTROL_HEADERS
        include/tdscontrol/tds.hpp
        include/tdscontrol/roots.hpp
        src/newtons_method.hpp
        src/discretization/spectral.hpp
        src/compute_N_rhp.hpp
        src/utils/ipow.hpp
        src/utils/thomas_algorithm.hpp
        src/utils/cubic_spline.hpp
        src/utils/grid_iterator.hpp
        tests/convert_reader_data.hpp)

set(TDS_CONTROL_SRC src/tds.cpp src/roots.cpp src/discretization/spectral.cpp
        src/compute_N_rhp.cpp
        src/utils/grid_iterator.cpp
        src/utils/cubic_spline.cpp
        src/utils/thomas_algorithm.cpp)

find_package(Eigen3 QUIET NO_MODULE)
if(NOT Eigen3_FOUND)
  message(STATUS "Eigen not found, fetching with FetchContent")

  set(EIGEN_BUILD_TESTS
      OFF
      CACHE BOOL "" FORCE)
  FetchContent_Declare(
    eigen
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    GIT_TAG 3.4.0)

  FetchContent_MakeAvailable(eigen)

endif()

add_library(tdscontrol STATIC ${TDS_CONTROL_HEADERS} ${TDS_CONTROL_SRC})
set_property(TARGET tdscontrol PROPERTY POSITION_INDEPENDENT_CODE ON)
target_include_directories(tdscontrol PUBLIC include)
target_link_libraries(tdscontrol PUBLIC Eigen3::Eigen)
target_compile_features(tdscontrol PRIVATE cxx_std_20)

if(BUILD_EXAMPLES)
  add_executable(demo examples/demo.cpp)
  target_link_libraries(demo tdscontrol)
endif()

if(TDS_BUILD_TESTING)
  message(STATUS "Build tests")
  include(CTest)
  include(Sanitizers)
  include(Warnings)
  include(AddTest)

  enable_testing()
  find_package(GTest REQUIRED)

  add_subdirectory(TDSControl-examples)

  enable_sanitizers(tdscontrol)
  enable_warnings(tdscontrol)
  target_compile_definitions(tdscontrol PUBLIC TDS_TESTING)

  add_tdscontrol_test(test_roots "tests/test_roots.cpp")
  add_tdscontrol_test(test_newton "tests/newtons_method.cpp")
  add_tdscontrol_test(test_sparse_spectral_discretization
                      "tests/sparse_spectral_discretization.cpp")
  add_tdscontrol_test(test_grid_iterator "tests/utils/grid_iterator.cpp")
  add_tdscontrol_test(test_splines "tests/utils/cubic_spline.cpp")
  add_tdscontrol_test(test_compute_N_rhp "tests/compute_N_rhp.cpp")
endif()

find_program(CLANG_FORMAT clang-format)
if(CLANG_FORMAT)
  # TODO: make sure that this also catches new files!
  file(
    GLOB_RECURSE
    ALL_CXX_SOURCE_FILES
    CONFIGURE_DEPENDS
    src/*.[ch]pp
    include/*.hpp
    tests/*.[ch]pp
    examples/*.[ch]pp)

  add_custom_target(
    format
    COMMAND ${CLANG_FORMAT} -i ${ALL_CXX_SOURCE_FILES}
    COMMENT "Running clang-format on source files")
endif()

set(CMAKE_CXX_CLANG_TIDY clang-tidy)
if(CMAKE_CXX_CLANG_TIDY)
  add_custom_target(
    tidy
    COMMAND run-clang-tidy -p ${CMAKE_BINARY_DIR}
    COMMENT "Running clang-tidy on the entire project")
endif()
