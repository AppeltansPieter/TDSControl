function(enable_warnings target)
  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")

    target_compile_options(${target} PRIVATE -Wall -Wextra -Wpedantic)

    if(TDS_ENABLE_WARNINGS_AS_ERRORS)
      target_compile_options(${target} PRIVATE -Werror)
    endif()

  elseif(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")

    target_compile_options(${target} PRIVATE /W4 /permissive-)

    if(TDS_ENABLE_WARNINGS_AS_ERRORS)
      target_compile_options(${target} PRIVATE /WX)
    endif()

  endif()
endfunction()
