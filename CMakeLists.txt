cmake_minimum_required(VERSION 3.16.0)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
option(TDS_ENABLE_COVERAGE "Enable coverage reporting" OFF)
option(TDS_ENABLE_WARNINGS_AS_ERRORS "Treat warnings as errors" OFF)
option(TDS_BUILD_TESTING "Build the tests" OFF)
option(TDS_BUILD_EXAMPLES "Build the examples" OFF)

file(READ "${CMAKE_CURRENT_SOURCE_DIR}/VERSION" VERSION_CONTENTS)
string(REGEX MATCH "VERSION = \"([0-9]+\.[0-9]+\.[0-9])\"" _
             ${VERSION_CONTENTS})
set(PROJECT_VERSION "${CMAKE_MATCH_1}")

project(
  tds-control
  VERSION ${PROJECT_VERSION}
  LANGUAGES C CXX)

if(TDS_ENABLE_COVERAGE)
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        message(STATUS "Coverage enabled")
        add_compile_options(-O0 -g --coverage)
        add_link_options(--coverage)
    endif()
endif()

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
target_include_directories(tdscontrol PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>)
target_link_libraries(tdscontrol PUBLIC Eigen3::Eigen)
target_compile_features(tdscontrol PRIVATE cxx_std_20)
add_library(TDSControl::tdscontrol ALIAS tdscontrol)

if(TDS_BUILD_EXAMPLES)
  add_subdirectory(examples)
endif()

if(TDS_BUILD_TESTING)
  add_subdirectory(TDSControl-examples)
  add_subdirectory(tests)
endif()

## Installation
include(GNUInstallDirs)
### Install the tdscontrol static library
install(TARGETS tdscontrol
        EXPORT TDSControlTargets
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
### Install the header files
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

### Generate TDSControlTargets.cmake file and install it.
install(
  EXPORT TDSControlTargets
  FILE
    TDSControlTargets.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/TDSControl
  NAMESPACE TDSControl::
)

### Create TDSControlConfig.cmake and TDSControlVersion.cmake and install these
include(CMakePackageConfigHelpers)
configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/TDSControlConfig.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/TDSControlConfig.cmake"
  INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/TDSControl"
  NO_SET_AND_CHECK_MACRO
  NO_CHECK_REQUIRED_COMPONENTS_MACRO
  )

write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/TDSControlConfigVersion.cmake
        VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMinorVersion)
install(FILES
        ${CMAKE_CURRENT_BINARY_DIR}/TDSControlConfig.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/TDSControlConfigVersion.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/TDSControl
)

find_program(CLANG_FORMAT clang-format)
if(CLANG_FORMAT)
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

if(TDS_ENABLE_COVERAGE)
    find_program(GCOVR gcovr)
    if(GCOVR)
        add_custom_target(coverage_html
            COMMAND ${GCOVR} --root ${CMAKE_SOURCE_DIR}
                             --filter 'src/.*' --filter 'include/.*'
                             --html --html-details -o ${CMAKE_CURRENT_BINARY_DIR}/coverage.html
                             --object-directory ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles
                             --exclude-throw-branches
                             --exclude-unreachable-branches
                             --exclude-function-lines
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
            COMMENT "Generating coverage report..."
        )
        add_custom_target(coverage_json
            COMMAND ${GCOVR} -r ${CMAKE_SOURCE_DIR}
                             --filter 'src/.*' --filter 'include/.*'
                             --object-directory ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles
                             --exclude-throw-branches
                             --exclude-unreachable-branches
                             --exclude-function-lines
                             --json -o ${CMAKE_CURRENT_BINARY_DIR}/coverage.json
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
            COMMENT "Generating coverage report..."
        )
    else()
        message(STATUS "Gcovr not found!")
    endif()
endif()
