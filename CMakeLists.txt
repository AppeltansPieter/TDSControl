cmake_minimum_required(VERSION 3.16.0)
project(tds-control VERSION 0.0.1 LANGUAGES C CXX)
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

if (PYTHON_BINDINGS)
    set(TDS_CONTROL_PY_BINDINGS
            python/py_bindings.cpp)

    # Add Python extension build using nanobind
    find_package(Python 3.9 REQUIRED COMPONENTS Interpreter Development.Module REQUIRED)
    execute_process(
        COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
        OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE nanobind_ROOT)
    find_package(nanobind CONFIG REQUIRED)

    nanobind_add_module(pytdscontrol MODULE ${TDS_CONTROL_PY_BINDINGS})
    target_link_libraries(pytdscontrol PRIVATE tdscontrol)

    # Install the compiled extension into the package dir inside the wheel
    install(TARGETS pytdscontrol
            LIBRARY DESTINATION "tdscontrol"    # Unix .so
            RUNTIME DESTINATION "tdscontrol")   # Windows .pyd/.dll
endif()

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

if (OCTAVE_BINDINGS)
    find_program(MKOCTFILE_EXECUTABLE mkoctfile)

    if(NOT MKOCTFILE_EXECUTABLE)
        message(FATAL_ERROR "mkoctfile not found")
    endif()

    get_target_property(TDSCONTROL_INCLUDES
        tdscontrol
        INTERFACE_INCLUDE_DIRECTORIES
    )

    if(TDSCONTROL_INCLUDES)
        list(TRANSFORM TDSCONTROL_INCLUDES PREPEND "-I")
    endif()

    add_custom_command(
        OUTPUT tds_create.oct
        COMMAND ${MKOCTFILE_EXECUTABLE}
                -std=c++17
                ${TDSCONTROL_INCLUDES}
                ${CMAKE_CURRENT_SOURCE_DIR}/octave/tds_create.cpp
                $<TARGET_FILE:tdscontrol>
                -o tds_create.oct
        DEPENDS tdscontrol octave/wrap_tds_create.cpp
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )

    add_custom_command(
        OUTPUT tds_roots.oct
        COMMAND ${MKOCTFILE_EXECUTABLE}
                -std=c++17
                ${TDSCONTROL_INCLUDES}
                ${CMAKE_CURRENT_SOURCE_DIR}/octave/wrap_tds_roots.cpp
                $<TARGET_FILE:tdscontrol>
                -o tds_roots.oct
        DEPENDS tdscontrol ${CMAKE_CURRENT_SOURCE_DIR}/octave/tds_roots.cpp
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )

    add_custom_target(octave_bindings ALL
        DEPENDS tds_roots.oct tds_create.oct
    )
    #TODO: find octave location and install in packages!
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