# Compiler flags for every Spp target.
#
# Precision policy: speed must never cost accuracy.
#   - No fast-math (or any of its sub-flags): IEEE 754 semantics are required.
#   - No floating-point contraction: fusing a*b+c into FMA makes last-bit
#     results depend on the target CPU. Results must be reproducible across
#     machines, so contraction stays off.
#   - No hidden extended-precision intermediates.
#
# Developed and tested with GCC. The Clang and MSVC branches select the
# equivalent options but are untested.

function(spp_apply_flags target)
    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        target_compile_options(${target} PRIVATE
            -Wall -Wextra -Wpedantic -Wshadow -Wdouble-promotion
            -ffp-contract=off -fexcess-precision=standard -fno-fast-math)
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        target_compile_options(${target} PRIVATE
            -Wall -Wextra -Wpedantic -Wshadow -Wdouble-promotion
            -ffp-contract=off -fno-fast-math)
    elseif(MSVC)
        target_compile_options(${target} PRIVATE /W4 /permissive- /fp:precise /utf-8)
    endif()

    if(SPP_NATIVE AND NOT MSVC)
        target_compile_options(${target} PRIVATE -march=native)
    endif()
endfunction()
