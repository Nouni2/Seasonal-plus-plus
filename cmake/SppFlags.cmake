# Compiler flags for every Spp target.
#
# Precision policy: speed must never cost accuracy.
#   - No -ffast-math (or any of its sub-flags): IEEE 754 semantics are required.
#   - -ffp-contract=off: GCC in GNU mode fuses a*b+c into FMA when -march allows
#     it, which makes last-bit results depend on the target CPU. Results must
#     be reproducible across machines, so contraction stays off.
#   - -fexcess-precision=standard: no hidden extended-precision intermediates.

function(spp_apply_flags target)
    target_compile_options(${target} PRIVATE
        -Wall -Wextra -Wpedantic -Wshadow -Wdouble-promotion
        -ffp-contract=off
        -fexcess-precision=standard
        -fno-fast-math
    )
    if(SPP_NATIVE)
        target_compile_options(${target} PRIVATE -march=native)
    endif()
endfunction()
