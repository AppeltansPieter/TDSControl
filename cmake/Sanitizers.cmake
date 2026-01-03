function(enable_sanitizers target)
    if (CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
        target_compile_options(${target} PRIVATE
            -fsanitize=address,undefined
            -fno-omit-frame-pointer
        )
        target_link_options(${target} PUBLIC
            -fsanitize=address,undefined
        )
    elseif (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
        message(WARNING "Sanitizers are not supported on MSVC. Skipping.")
    endif()
endfunction()