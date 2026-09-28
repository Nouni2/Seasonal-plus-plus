#pragma once

// Internal: build mode. With SPP_EXACT (reference mode), every series is
// evaluated directly and every formula keeps its textbook operation order,
// with the fixed-tolerance event solver. Otherwise the faster, equally
// accurate forms are used: Chebyshev caching, Horner polynomials, direct
// heliocentric-to-geocentric inversion, split GMST and the default solver.

namespace spp::detail {

#ifdef SPP_EXACT
inline constexpr bool kExact = true;
#else
inline constexpr bool kExact = false;
#endif

}  // namespace spp::detail
