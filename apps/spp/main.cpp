// spp: command-line front end for the Spp library.
//
//   spp sun  YYYY-MM-DD HH:MM[:SS] [LAT LON [ELEV]]   apparent Sun (UTC input)
//   spp day  YYYY-MM-DD LAT LON [ELEV]                transit, sunrise, sunset
//   spp year YYYY LAT LON [ELEV]                      one CSV row per day
//
// Options: --dut1=S (default: estimated)  --pressure=MBAR (default: standard at
// the elevation)  --temp=C  --no-refraction
// Longitude is positive East; elevation in metres.

#include <charconv>
#include <cstdio>
#include <cstdlib>
#include <optional>
#include <print>
#include <string>
#include <string_view>
#include <system_error>
#include <vector>

#include "spp/spp.hpp"

namespace {

struct Options {
    std::optional<double> dut1;   // estimated when not given
    spp::Atmosphere atm;
};

[[noreturn]] void usage() {
    std::println(stderr,
                 "usage:\n"
                 "  spp sun  YYYY-MM-DD HH:MM[:SS] [LAT LON [ELEV]]\n"
                 "  spp day  YYYY-MM-DD LAT LON [ELEV]\n"
                 "  spp year YYYY LAT LON [ELEV]\n"
                 "options: --dut1=S --pressure=MBAR --temp=C --no-refraction\n"
                 "  DUT1 is estimated and pressure follows the elevation unless given");
    std::exit(2);
}

double num(std::string_view s) {
    double v{};
    auto [p, ec] = std::from_chars(s.data(), s.data() + s.size(), v);
    if (ec != std::errc{} || p != s.data() + s.size()) {
        std::println(stderr, "spp: not a number: '{}'", s);
        std::exit(2);
    }
    return v;
}

int integer(std::string_view s) {
    int v{};
    auto [p, ec] = std::from_chars(s.data(), s.data() + s.size(), v);
    if (ec != std::errc{} || p != s.data() + s.size()) {
        std::println(stderr, "spp: not an integer: '{}'", s);
        std::exit(2);
    }
    return v;
}

// Splits "a<sep>b<sep>c" into integers/numbers; leading '-' allowed on the first field.
std::vector<std::string_view> fields(std::string_view s, char sep) {
    std::vector<std::string_view> out;
    std::size_t start = s.starts_with('-') ? 1 : 0;
    std::size_t from = 0;
    for (std::size_t i = start; i <= s.size(); ++i) {
        if (i == s.size() || s[i] == sep) {
            out.push_back(s.substr(from, i - from));
            from = i + 1;
        }
    }
    return out;
}

spp::CivilDate parse_date(std::string_view s) {
    auto f = fields(s, '-');
    if (f.size() != 3) usage();
    return {integer(f[0]), integer(f[1]), integer(f[2])};
}

spp::Observer parse_observer(const std::vector<std::string_view>& a, std::size_t i) {
    return {num(a[i]), num(a[i + 1]), a.size() > i + 2 ? num(a[i + 2]) : 0.0};
}

// Full ISO timestamp: events can fall on the neighbouring UTC date.
std::string when(const std::optional<spp::Time>& t) {
    if (!t) return "-";
    return t->iso_utc();
}

const char* type_name(spp::events::DayType t) {
    switch (t) {
        case spp::events::DayType::Normal: return "NORMAL";
        case spp::events::DayType::PolarDay: return "POLAR_DAY";
        case spp::events::DayType::PolarNight: return "POLAR_NIGHT";
    }
    return "?";
}

int cmd_sun(const std::vector<std::string_view>& a, const Options& o) {
    if (a.size() != 2 && a.size() != 4 && a.size() != 5) usage();
    const auto d = parse_date(a[0]);
    const auto hms = fields(a[1], ':');
    if (hms.size() < 2 || hms.size() > 3) usage();
    const auto t = spp::Time::from_utc(d.year, d.month, d.day, integer(hms[0]), integer(hms[1]),
                                       hms.size() == 3 ? num(hms[2]) : 0.0, o.dut1);

    const auto tt = t.tt();
    const auto n = spp::earth::nutation(tt);
    const auto ecl = spp::sun::apparent_ecliptic(tt, n);
    const auto eq = spp::sun::apparent(tt, n);

    std::println("time        {}  (ΔT {:.2f} s, DUT1 {:+.3f} s{})", t.iso_utc(), t.delta_t(), t.dut1(),
                 t.dut1_estimated() ? " estimated" : "");
    std::println("ecliptic    λ {:.8f}°  β {:.8f}°  R {:.9f} AU", ecl.lon_deg, ecl.lat_deg, ecl.dist_au);
    std::println("equatorial  α {:.8f}° ({:.8f} h)  δ {:.8f}°", eq.ra_deg, eq.ra_hours(), eq.dec_deg);
    std::println("nutation    Δψ {:.4f}\"  Δε {:.4f}\"  ε {:.8f}°", n.dpsi_deg * 3600.0,
                 n.deps_deg * 3600.0, n.true_obliquity_deg);
    std::println("sidereal    GMST {:.8f}°  GAST {:.8f}°", spp::earth::gmst_deg(t.ut1()),
                 spp::earth::gast_deg(t.ut1(), n));

    if (a.size() >= 4) {
        const auto obs = parse_observer(a, 2);
        const auto hor = spp::sun::horizontal(t, obs, o.atm);
        const auto geo = spp::sun::horizontal(t, obs, {.refraction = false});
        std::println("horizontal  alt {:.6f}° (geometric {:.6f}°)  az {:.6f}°", hor.alt_deg, geo.alt_deg, hor.az_deg);
    }
    return 0;
}

int cmd_day(const std::vector<std::string_view>& a, const Options& o) {
    if (a.size() != 3 && a.size() != 4) usage();
    const auto d = spp::events::day(parse_date(a[0]), parse_observer(a, 1), o.atm, o.dut1);
    std::println("type      {}", type_name(d.type));
    std::println("transit   {}", d.transit.iso_utc());
    std::println("sunrise   {}", d.rise ? d.rise->iso_utc() : "-");
    std::println("sunset    {}", d.set ? d.set->iso_utc() : "-");
    std::println("duration  {:.4f} h", d.duration_hours);
    return 0;
}

bool leap(int y, bool gregorian) {
    return gregorian ? (y % 4 == 0 && (y % 100 != 0 || y % 400 == 0)) : (y % 4 == 0);
}

int cmd_year(const std::vector<std::string_view>& a, const Options& o) {
    if (a.size() != 3 && a.size() != 4) usage();
    const int year = integer(a[0]);
    const auto obs = parse_observer(a, 1);
    static constexpr int kDays[] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};

    std::println("date,type,transit_utc,sunrise_utc,sunset_utc,duration_h");
    for (int m = 1; m <= 12; ++m) {
        const int n = kDays[m - 1] + (m == 2 && leap(year, year > 1582) ? 1 : 0);
        for (int dd = 1; dd <= n; ++dd) {
            const spp::CivilDate date{year, m, dd};
            if (year == 1582 && m == 10 && dd > 4 && dd < 15) continue;   // dropped by the reform
            const auto d = spp::events::day(date, obs, o.atm, o.dut1);
            std::println("{:04d}-{:02d}-{:02d},{},{},{},{},{:.6f}", year, m, dd, type_name(d.type),
                         when(d.transit), when(d.rise), when(d.set), d.duration_hours);
        }
    }
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    Options o;
    std::vector<std::string_view> args;
    for (int i = 1; i < argc; ++i) {
        std::string_view s = argv[i];
        if (s.starts_with("--dut1=")) o.dut1 = num(s.substr(7));
        else if (s.starts_with("--pressure=")) o.atm.pressure_mbar = num(s.substr(11));
        else if (s.starts_with("--temp=")) o.atm.temp_c = num(s.substr(7));
        else if (s == "--no-refraction") o.atm.refraction = false;
        else args.push_back(s);
    }
    if (args.empty()) usage();

    const std::string_view cmd = args.front();
    args.erase(args.begin());
    try {
        if (cmd == "sun") return cmd_sun(args, o);
        if (cmd == "day") return cmd_day(args, o);
        if (cmd == "year") return cmd_year(args, o);
    } catch (const std::system_error&) {
        return 1;   // output closed early, e.g. piped into `head`
    }
    usage();
}
