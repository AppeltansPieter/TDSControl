function(add_tdscontrol_test executable_name file_names)
  add_executable("${executable_name}" "${file_names}")
  target_link_libraries(${executable_name} PRIVATE GTest::gtest_main
                                                   GTest::gmock tdscontrol tdsreader)
  enable_sanitizers(${executable_name})
  enable_warnings(${executable_name})
  gtest_discover_tests(${executable_name}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
          PROPERTIES TIMEOUT 10 ENVIRONMENT "GTEST_COLOR=1"
  )
endfunction()
