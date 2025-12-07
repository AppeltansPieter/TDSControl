cmake_minimum_required(VERSION 3.16.0)
project(tds-control VERSION 0.0.1 LANGUAGES C CXX)

include(CTest)

set(TDS_CONTROL_HEADERS
        include/tdscontrol/tds.hpp
        include/tdscontrol/roots.hpp)

set(TDS_CONTROL_SRC 
            src/tds.cpp
            src/roots.cpp)

set(TDS_CONTROL_TEST
            tests/test_roots.cpp)

find_package(Eigen3 REQUIRED NO_MODULE)

add_library(tdscontrol SHARED ${TDS_CONTROL_HEADERS} ${TDS_CONTROL_SRC})
target_include_directories(tdscontrol PUBLIC include)
target_link_libraries(tdscontrol PRIVATE Eigen3::Eigen)

add_executable(demo examples/demo.cpp)
target_link_libraries(demo tdscontrol)

if (BUILD_TESTING)
    message(STATUS "Build tests")
    enable_testing()
    find_package(GTest REQUIRED)

    add_executable(test_tds ${TDS_CONTROL_TEST})
    target_link_libraries(test_tds PRIVATE GTest::gtest_main GTest::gmock tdscontrol)
    add_test(NAME test_tds COMMAND test_tds)
endif()