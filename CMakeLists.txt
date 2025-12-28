cmake_minimum_required(VERSION 3.16.0)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
option(ENABLE_WARNINGS_AS_ERRORS "Treat warnings as errors" OFF)

file(READ "${CMAKE_SOURCE_DIR}/VERSION" VERSION_CONTENTS)
string(REGEX MATCH "VERSION = \"([0-9]+\.[0-9]+\.[0-9])\"" _
             ${VERSION_CONTENTS})
set(PROJECT_VERSION "${CMAKE_MATCH_1}")

project(
  tds-control
  VERSION ${PROJECT_VERSION}
  LANGUAGES C CXX)

include(FetchContent)

set(TDS_CONTROL_HEADERS include/tdscontrol/tds.hpp include/tdscontrol/roots.hpp)

set(TDS_CONTROL_SRC src/tds.cpp src/roots.cpp)

set(TDS_CONTROL_TEST tests/test_roots.cpp)

find_package(Eigen3 QUIET NO_MODULE)
if(NOT Eigen3_FOUND)
  message(STATUS "Eigen not found, fetching with FetchContent")

  FetchContent_Declare(
    eigen
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    GIT_TAG 3.4.0)

  FetchContent_MakeAvailable(eigen)

endif()

add_library(tdscontrol STATIC ${TDS_CONTROL_HEADERS} ${TDS_CONTROL_SRC})
set_property(TARGET tdscontrol PROPERTY POSITION_INDEPENDENT_CODE ON)
target_include_directories(tdscontrol PUBLIC include)
target_link_libraries(tdscontrol PRIVATE Eigen3::Eigen)

if(BUILD_EXAMPLES)
  add_executable(demo examples/demo.cpp)
  target_link_libraries(demo tdscontrol)
endif()

if(BUILD_TESTING)
  message(STATUS "Build tests")
  include(CTest)
  include(Sanitizers)
  include(Warnings)

  enable_testing()
  find_package(GTest REQUIRED)

  add_executable(test_tds ${TDS_CONTROL_TEST})
  target_link_libraries(test_tds PRIVATE GTest::gtest_main GTest::gmock
                                         tdscontrol)
  enable_sanitizers(test_tds)
  enable_sanitizers(tdscontrol)

  enable_warnings(test_tds)
  enable_warnings(tdscontrol)
  add_test(NAME test_tds COMMAND test_tds)
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
