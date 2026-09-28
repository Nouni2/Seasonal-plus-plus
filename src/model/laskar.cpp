#include "spp/model/laskar.hpp"

#include <cmath>

#include "../build_mode.hpp"

namespace spp::model::laskar {

double mean_obliquity_deg(double T) {
    const double U = T / 100.0;
    if constexpr (!detail::kExact) {
        // Horner form of the same polynomial (arcseconds).
        const double p = 2.45;
        const double q = ((((((((p * U + 5.79) * U + 27.87) * U + 7.12) * U - 39.05) * U - 249.67) * U - 51.38)
                            * U + 1999.25) * U - 1.55) * U - 4680.93;
        return (84381.448 + q * U) / 3600.0;
    }
    const double eps0_arcsec =
        84381.448
        - 4680.93 * U
        - 1.55 * std::pow(U, 2.0)
        + 1999.25 * std::pow(U, 3.0)
        - 51.38 * std::pow(U, 4.0)
        - 249.67 * std::pow(U, 5.0)
        - 39.05 * std::pow(U, 6.0)
        + 7.12 * std::pow(U, 7.0)
        + 27.87 * std::pow(U, 8.0)
        + 5.79 * std::pow(U, 9.0)
        + 2.45 * std::pow(U, 10.0);
    return eps0_arcsec / 3600.0;
}

}  // namespace spp::model::laskar
