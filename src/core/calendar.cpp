#include "spp/core/calendar.hpp"

namespace spp {

bool is_gregorian(const CivilDate& d) {
    if (d.year != 1582) return d.year > 1582;
    if (d.month != 10) return d.month > 10;
    return d.day >= 15;
}

}  // namespace spp
