cmake_minimum_required(VERSION 3.16.0)

file(READ "${CMAKE_SOURCE_DIR}/VERSION" VERSION_CONTENTS) 
string(REGEX MATCH "VERSION = \"([0-9]+\.[0-9]+\.[0-9])\"" _ ${VERSION_CONTENTS})
set(PROJECT_VERSION "${CMAKE_MATCH_1}")

project(tds-control VERSION ${PROJECT_VERSION} LANGUAGES C CXX)
include(FetchContent)

set(TDS_CONTROL_HEADERS
        include/tdscontrol/tds.hpp
        include/tdscontrol/roots.hpp)

set(TDS_CONTROL_SRC 
            src/tds.cpp
            src/roots.cpp)

set(TDS_CONTROL_TEST
            tests/test_roots.cpp)

find_package(Eigen3 QUIET NO_MODULE)
if(NOT Eigen3_FOUND)
  message(STATUS "Eigen not found, fetching with FetchContent")

  FetchContent_Declare(
    eigen
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    GIT_TAG 3.4.0   
  )

  FetchContent_MakeAvailable(eigen)

endif()

add_library(tdscontrol STATIC ${TDS_CONTROL_HEADERS} ${TDS_CONTROL_SRC})
set_property(TARGET tdscontrol PROPERTY POSITION_INDEPENDENT_CODE ON)
target_include_directories(tdscontrol PUBLIC include)
target_link_libraries(tdscontrol PRIVATE Eigen3::Eigen)

add_executable(demo examples/demo.cpp)
target_link_libraries(demo tdscontrol)


if (JULIA_BINDINGS)
  find_package(JlCxx)
  get_target_property(JlCxx_location JlCxx::cxxwrap_julia LOCATION)
  get_filename_component(JlCxx_location ${JlCxx_location} DIRECTORY)
  set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib;${JlCxx_location}")

  message(STATUS "Found JlCxx at ${JlCxx_location}")
  set_target_properties(JlCxx::cxxwrap_julia PROPERTIES INTERFACE_COMPILE_FEATURES "")
  set(CMAKE_CXX_STANDARD 17)

  target_compile_definitions(tdscontrol PUBLIC JULIA_BINDING)
  target_link_libraries(tdscontrol PUBLIC JlCxx::cxxwrap_julia)

  add_library(tdscontroljl SHARED julia/bindings.cpp)
  target_link_libraries(tdscontroljl PUBLIC JlCxx::cxxwrap_julia tdscontrol)

  install(TARGETS tdscontroljl 
          LIBRARY DESTINATION lib
          ARCHIVE DESTINATION lib
          RUNTIME DESTINATION lib)
endif()

if (BUILD_TESTING)
    message(STATUS "Build tests")
    include(CTest)
    enable_testing()
    find_package(GTest REQUIRED)

    add_executable(test_tds ${TDS_CONTROL_TEST})
    target_link_libraries(test_tds PRIVATE GTest::gtest_main GTest::gmock tdscontrol)
    add_test(NAME test_tds COMMAND test_tds)
endif()