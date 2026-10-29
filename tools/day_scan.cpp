// day_scan: consistency of events::day with a dense sampling of its own condition.
//
// For every day of a year at |latitude| from 60° to 90° in 0.25° steps, both
// hemispheres, the apparent altitude of the upper limb is sampled every 2
// minutes over the 24 h centered on transit and compared with the result of
// events::day:
//   - PolarNight: the limb is below the horizon at every sample;
//   - PolarDay:   the limb is above the horizon at every sample;
//   - Normal:     at least one of sunrise and sunset is present, each where the
//                 limb's altitude is 0 within 1e-5", and the samples change sign.
// Prints the inconsistencies, PASS or FAIL, and exits with 0 or 1.
//
//   spp_day_scan [YEAR]   (default 2026)

#include <charconv>
#include <cmath>
#include <print>
#include <string_view>

#include "spp/spp.hpp"

int main(int argc, char** argv) {
    int year = 2026;
    if (argc > 1) std::from_chars(argv[1], argv[1] + std::string_view(argv[1]).size(), year);
    long days = 0, bad = 0;
    for (int hemi : {1, -1}) {
        for (int k = 0; k <= 120; ++k) {
            const spp::Observer obs{hemi * (60.0 + 0.25 * k), 0.0, 0.0};
            auto t = spp::Time::from_utc(year, 1, 1, 12);
            while (t.to_civil().year == year) {
                const auto c = t.to_civil();
                const spp::CivilDate date{c.year, c.month, c.day};
                const auto d = spp::events::day(date, obs);
                int up = 0, down = 0;
                for (int s = -360; s <= 360; ++s) {
                    const double u = spp::events::upper_limb_altitude_deg(d.transit.shifted(s / 720.0), obs, {});
                    (u > 0.0 ? up : down) += 1;
                }
                const char* why = nullptr;
                using spp::events::DayType;
                if (d.type == DayType::PolarNight && up > 0) why = "PolarNight, but the limb is up at some sample";
                if (d.type == DayType::PolarDay && down > 0) why = "PolarDay, but the limb is down at some sample";
                if (d.type == DayType::Normal) {
                    if (!d.rise && !d.set) why = "Normal without sunrise or sunset";
                    else if (up == 0 || down == 0) why = "Normal, but the samples never change sign";
                    for (const auto& e : {d.rise, d.set})
                        if (e && std::fabs(spp::events::upper_limb_altitude_deg(*e, obs, {})) * 3600.0 > 1e-5)
                            why = "event where the limb is not at 0";
                }
                if (why) {
                    ++bad;
                    if (bad <= 20) std::println("lat {:+.2f} {}: {}", obs.lat_deg, spp::iso_date(date), why);
                }
                ++days;
                t = t.shifted(1.0);
            }
        }
    }
    std::println("{} days, {} inconsistent", days, bad);
    std::println("{}", bad ? "FAIL" : "PASS");
    return bad ? 1 : 0;
}
