function(add_tdscontrol_test executable_name file_names)
  add_executable("${executable_name}" "${file_names}")
  target_link_libraries(${executable_name} PRIVATE GTest::gtest_main
                                                   GTest::gmock tdscontrol)
  enable_sanitizers(${executable_name})
  enable_warnings(${executable_name})
  add_test(NAME ${executable_name} COMMAND ${executable_name})
endfunction()
