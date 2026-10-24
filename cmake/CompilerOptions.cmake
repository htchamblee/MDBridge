# Compiler settings shared by every MDBridge target.

if(MSVC)
    # sources and string literals are UTF-8; __cplusplus reports the real standard; headers included
    # with <...> (Qt, the fetched libraries) are third-party code and give no warnings
    add_compile_options(/utf-8 /Zc:__cplusplus /permissive- /external:anglebrackets /external:W0)
    add_compile_definitions(NOMINMAX WIN32_LEAN_AND_MEAN _CRT_SECURE_NO_WARNINGS)
endif()

# Warnings for MDBridge's own targets (not for fetched dependencies).
function(mdbridge_target_warnings target)
    if(MSVC)
        target_compile_options(${target} PRIVATE /W4)
        if(MDBRIDGE_WARNINGS_AS_ERRORS)
            target_compile_options(${target} PRIVATE /WX)
        endif()
    else()
        target_compile_options(${target} PRIVATE -Wall -Wextra -Wpedantic -Wnon-virtual-dtor
                                                  -Wold-style-cast -Wcast-align -Woverloaded-virtual
                                                  -Wno-missing-field-initializers)  # Vulkan structs: omitted = zero
        if(MDBRIDGE_WARNINGS_AS_ERRORS)
            target_compile_options(${target} PRIVATE -Werror)
        endif()
    endif()
endfunction()
