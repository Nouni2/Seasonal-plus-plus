#include "spp/model/laskar.hpp"

namespace spp::model::laskar {

double mean_obliquity_deg(double T) {
    // ε0 = 84381.448" - 4680.93"U - 1.55"U² + 1999.25"U³ - 51.38"U⁴ - 249.67"U⁵
    //      - 39.05"U⁶ + 7.12"U⁷ + 27.87"U⁸ + 5.79"U⁹ + 2.45"U¹⁰, in Horner form.
    const double U = T / 100.0;
    const double q = ((((((((2.45 * U + 5.79) * U + 27.87) * U + 7.12) * U - 39.05) * U - 249.67) * U - 51.38)
                        * U + 1999.25) * U - 1.55) * U - 4680.93;
    return (84381.448 + q * U) / 3600.0;
}

}  // namespace spp::model::laskar
