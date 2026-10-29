# Compiler flags for every Spp target.
#
# Precision policy: speed must never cost accuracy.
#   - No fast-math (or any of its sub-flags): IEEE 754 semantics are required.
#   - No floating-point contraction: the compiler never fuses a*b+c into an
#     FMA, so the arithmetic it generates does not depend on the target CPU.
#   - No hidden extended-precision intermediates. This matters only where
#     floating point runs on the x87 unit (32-bit x86 without SSE2); on
#     x86-64, which uses SSE2, the flag does not change the generated code.
#
# These make the compiler's code generation deterministic, not the results.
# sin, cos, atan2, pow and the other transcendental functions come from the
# platform's C math library (glibc, the Windows UCRT, Apple libm), which is not
# correctly rounded: results can differ in the last bits between platforms,
# and glibc on x86-64 also selects FMA/AVX2 variants by CPU at run time
# (README, Floating-point policy).
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
